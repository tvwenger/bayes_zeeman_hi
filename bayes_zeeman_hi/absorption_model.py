"""
absorption_model.py
AbsorptionModel definition

Copyright(C) 2025 by
Trey V. Wenger; tvwenger@gmail.com
This code is licensed under MIT license (see LICENSE for details)
"""

from typing import Iterable, Optional

import pymc as pm
import pytensor.tensor as pt

from bayes_spec import BaseModel

from bayes_zeeman_hi import physics

from typing import Iterable

import astropy.constants as c
import astropy.units as u


class AbsorptionModel(BaseModel):
    """Definition of the AbsorptionModel model."""

    def __init__(self, Tbg, *args, **kwargs):
        """Initialize a new AbsorptionModel instance,.

        Parameters
        ----------
        Tbg : float
            Background source brightness temperature (K)
        """
        # Initialize BaseModel
        super().__init__(*args, **kwargs)

        # Select features used for posterior clustering
        self._cluster_features += [
            "velocity",
            "fwhm2",
        ]

        # Define TeX representation of each parameter
        self.var_name_map.update(
            {
                "fwhm2": r"$\Delta V^2$ (km$^{2}$ s$^{-2}$)",
                "velocity": r"$V_{\rm LSR}$ (km s$^{-1}$)",
                "Bparallel": r"$B_{\rm los}$",
                "tau_total": r"$\int \tau(v) dv$ (km s$^{-1}$)",
            }
        )

        # Save inputs
        self.Tbg = Tbg

    def add_priors(
        self,
        prior_tau_total: float = 1.0,
        prior_fwhm2: float = 200.0,
        prior_velocity: Iterable[float] = [-10.0, 10.0],
        prior_Bparallel: Iterable[float] = [-20.0, 20.0],
        prior_leakage_fraction: float = 0.01,
    ):
        """Add priors and deterministics to the model

        Parameters
        ----------
        prior_tau_total : float, optional
            Prior distribution on total optical depth (km s-1), by default 1.0, where
            tau_total ~ HalfNormal(sigma=prior)
        prior_fwhm2 : float, optional
            Prior distribution on FWHM^2  (km2 s-2), by default 200.0, where
            fwhm2 ~ prior * ChiSquared(nu=1)
            i.e., half-normal on FWHM
        prior_velocity : Iterable[float], optional
            Prior distribution on centroid velocity (km s-1), by default [-10.0, 10.0], where
            velocity_norm ~ Beta(alpha=2.0, beta=2.0)
            velocity ~ prior[0] + (prior[1] - prior[0]) * velocity_norm
        prior_Bparallel : Iterable[float], optional
            Prior distribution on parallel magnetic field strength (uG), by default [-20.0, 20.0], where
            Bparallel_norm ~ Beta(alpha=2.0, beta=2.0)
            Bparallel ~ prior[0] + (prior[1] - prior[0]) * Bparallel_norm
        prior_leakage_fraction : float, optional
            Prior distribution on leakage fraction, by default 0.01, where
            leakage_fraction ~ HalfNormal(sigma=prior)
        """
        with self.model:
            # total optical depth (km s-1; shape: clouds)
            tau_total_norm = pm.HalfNormal("tau_total_norm", sigma=1.0, dims="cloud")
            _ = pm.Deterministic(
                "tau_total", prior_tau_total * tau_total_norm, dims="cloud"
            )

            # FWHM^2 (km2 s-2; shape: clouds)
            fwhm2_norm = pm.ChiSquared("fwhm2_norm", nu=1, dims="cloud")
            _ = pm.Deterministic("fwhm2", prior_fwhm2 * fwhm2_norm, dims="cloud")

            # Velocity (km/s; shape: clouds)
            velocity_norm = pm.Beta("velocity_norm", alpha=2.0, beta=2.0, dims="cloud")
            _ = pm.Deterministic(
                "velocity",
                prior_velocity[0]
                + (prior_velocity[1] - prior_velocity[0]) * velocity_norm,
                dims="cloud",
            )

            # Bparallel (uG; shape: clouds)
            Bparallel_norm = pm.Beta(
                "Bparallel_norm", alpha=2.0, beta=2.0, dims="cloud"
            )
            _ = pm.Deterministic(
                "Bparallel",
                prior_Bparallel[0]
                + (prior_Bparallel[1] - prior_Bparallel[0]) * Bparallel_norm,
                dims="cloud",
            )

            # leakage fraction
            leakage_fraction_norm = pm.HalfNormal("leakage_fraction_norm", sigma=1.0)
            _ = pm.Deterministic(
                "leakage_fraction", prior_leakage_fraction * leakage_fraction_norm
            )

    def add_likelihood(self):
        """Add likelihood to the model. SpecData key must be "I" and "V"."""
        # Evaluate line profile (shape: spectral, clouds)
        line_profile = physics.gaussian(
            self.data["I"].spectral[:, None],
            self.model["velocity"],
            pt.sqrt(self.model["fwhm2"]),
        )

        # LCP, RCP line profile (shape: spectra, clouds)
        z = (2.8 * u.Hz * c.c / (1420.4058 * u.MHz)).to("km/s").value  # km/s / uG
        line_profile_LCP = physics.gaussian(
            self.data["V"].spectral[:, None],
            self.model["velocity"] + z * self.model["Bparallel"] / 2.0,
            pt.sqrt(self.model["fwhm2"]),
        )
        line_profile_RCP = physics.gaussian(
            self.data["V"].spectral[:, None],
            self.model["velocity"] - z * self.model["Bparallel"] / 2.0,
            pt.sqrt(self.model["fwhm2"]),
        )

        # Optical depth spectra (shape: spectral, clouds)
        optical_depth = self.model["tau_total"] * line_profile
        optical_depth_LCP = self.model["tau_total"] * line_profile_LCP
        optical_depth_RCP = self.model["tau_total"] * line_profile_RCP

        # Sum over clouds
        stokesI = 2.0 * self.Tbg * pt.exp(-optical_depth.sum(axis=1))
        stokesV = (
            self.Tbg
            * (
                pt.exp(-optical_depth_RCP.sum(axis=1))
                - pt.exp(-optical_depth_LCP.sum(axis=1))
            )
            + self.model["leakage_fraction"] * stokesI
        )

        with self.model:
            # Evaluate likelihood
            _ = pm.Normal(
                "I",
                mu=stokesI,
                sigma=self.data["I"].noise,
                observed=self.data["I"].brightness,
            )
            _ = pm.Normal(
                "V",
                mu=stokesV,
                sigma=self.data["V"].noise,
                observed=self.data["V"].brightness,
            )
