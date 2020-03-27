# General python imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
import yaml
import pickle
from pathlib import Path
import scipy.integrate as integrate
# GRB population specific imports
import stats as st
import physics as ph
import functional_forms as ff
import f90_functions as f90f
from io_grb_pop import root_dir, load_observational_constraints
from cosmology import init_cosmology, Lum_dist
from constants import R_tot_BATSE, R_tot_BATSE_err, T_live_BATSE

log = logging.getLogger(__name__)


def create_GRB_population_from(Nb_GRBs, cosmo, params, incl_instruments, incl_samples, incl_constraints,
    ECLAIRs_prop=None, output_dir=None, nGRB0=None, nGRB0_err=None, run_mode=None, savefig=False):
    """
        A convenienve function to quickly create a GRB population given
        a set of cosmology, parameters, instruments, samples,
        constraints and a number of GRBs to draw.
    """

    GRB_pop = GRBPopulation(Nb_GRBs=Nb_GRBs,
                            output_dir=output_dir)
    GRB_pop.draw_GRB_properties(cosmo=cosmo, params=params, run_mode=run_mode, savefig=savefig)
    GRB_pop.calculate_quantities(instruments=incl_instruments, samples=incl_samples, ECLAIRs_prop=ECLAIRs_prop)
    GRB_pop.create_mock_constraints(constraints=incl_constraints)
    GRB_pop.compare_to_observational_constraints(constraints=incl_constraints)
    if nGRB0 is None:
        GRB_pop.normalize_to_Stern()
    else:
        GRB_pop.normalize_from(nGRB0)
    print(GRB_pop.summary())
    return GRB_pop


class GRBPopulation:
    def __init__(self, Nb_GRBs=None, properties=None, output_dir=None):

        if Nb_GRBs is None and properties is None:
            raise ValueError("You must either specify Nb_GRBs or properties to instanciate this class.")
        elif Nb_GRBs is not None and properties is None:
            if isinstance(Nb_GRBs, str):
                Nb_GRBs = int(float(Nb_GRBs))
            if not isinstance(Nb_GRBs, int):
                raise TypeError
            self.Nb_GRBs = Nb_GRBs
            self.properties = pd.DataFrame({})
        elif Nb_GRBs is None and properties is not None:
            if not isinstance(properties, pd.DataFrame):
                raise ValueError("Properties must be a pandas DataFrame instance.")
            self.Nb_GRBs = len(properties)
            self.properties = properties
        self.parameters = {}
        self.mock_constraints = {}
        self.likelihood_params = {}
        self.normalization = {}
        if output_dir is None:
            self.output_dir = Path().cwd()
        else:
            if not isinstance(output_dir, Path):
                output_dir = Path(output_dir)
            self.output_dir = output_dir

    def draw_L(self, Nb_GRBs=None, model='EPL', z=None, run_mode=None, savefig=False, **params):
        """
            Draw L from the desired luminosity function
        """

        self.parameters['luminosity_function'] = {'model':model, **params}

        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if model == 'EPL':
            logLmin = params['logLmin']
            logLmax = params['logLmax']
            slope = params['slope']
            k_evol = params['k_evol']
            t = np.random.rand(Nb_GRBs)

            L = 10.**logLmin * (1. - (1. - (10.**(logLmax-logLmin))**(1.-slope))*t)**(1. / (1.-slope))

            if k_evol != 0.0:
                if z is None:
                    if 'z' not in self.properties.keys():
                        raise KeyError('Could not find z to use in draw_L')
                    z = self.properties['z']
                L *= (1 + z)**k_evol

        elif model == 'ES':
            logLmin = params['logLmin']
            logLbreak = params['logLbreak']
            logLmax = logLbreak + 1.5  # to avoid extremely low probabilities beyond logLbreak + 1.5
            slope = params['slope']
            k_evol = params['k_evol']
            logL_range = np.linspace(logLmin, logLmax, 1000)
            logL_pdf_unnormed = ff.Schechter_log(logL_range, logLbreak=logLbreak, slope=slope)
            logL_pdf = logL_pdf_unnormed / logL_pdf_unnormed.sum()

            L = 10.**(np.random.choice(logL_range, size=Nb_GRBs, p=logL_pdf))

            if k_evol != 0.0:
                if z is None:
                    if 'z' not in self.properties.keys():
                        raise KeyError('Could not find z to use in draw_L')
                    z = self.properties['z']
                L *= (1 + z)**k_evol

        elif model == 'EBPL':
            logLmin = params['logLmin']
            logLbreak = params['logLbreak']
            logLmax = params['logLmax']
            slopeL = params['slopeL']
            slopeH = params['slopeH']
            k_evol = params['k_evol']
            logL_range = np.linspace(logLmin, logLmax, 1000)
            logL_pdf_unnormed = ff.BPL_lum(logL_range, logLbreak=logLbreak, slopeL=slopeL, slopeH=slopeH)
            logL_pdf = logL_pdf_unnormed / logL_pdf_unnormed.sum()

            L = 10.**(np.random.choice(logL_range, size=Nb_GRBs, p=logL_pdf))

            if k_evol != 0.0:
                if z is None:
                    if 'z' not in self.properties.keys():
                        raise KeyError('Could not find z to use in draw_L')
                    z = self.properties['z']
                L *= (1 + z)**k_evol
        else:
            raise ValueError("Invalid model for L drawing")

        self.properties['L'] = L

        if run_mode == 'debug' and model != 'EPL':
            log.info("Debug mode activated; plotting L pdf")
            delta_logL = logL_range[1]-logL_range[0]
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            ax.hist(np.log10(L),
                    bins=np.arange(logLmin, logLmax, 0.1),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.plot(logL_range, logL_pdf/delta_logL, label='PDF')
            ax.set_xlabel('log(L [erg/s])')
            ax.set_ylabel('L PDF')
            ax.set_yscale('log')
            ax.legend()
            if savefig:
                self.save_fig(fig, 'L_draw.pdf')

        return L.copy()

    def draw_z(self, cosmo=None, Nb_GRBs=None, zmax=20, model='SH', run_mode=None, savefig=False, **params):
        """
            Draw z from a redshift distribution
        """

        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if cosmo is None:
            try:
                cosmo = init_cosmology(Path().cwd().parent/'data'/'cosmology')
            except FileNotFoundError:
                raise ValueError("Could not draw z because no cosmology provided")

        self.parameters['cosmology'] = {**cosmo['parameters']}
        self.parameters['redshift_distribution'] = {'model':model, **params}

        redshift = cosmo['redshift']
        dVdz = cosmo['dVdz']

        if zmax > redshift[-1]:
            log.warning("In draw_z, zmax is greater than largest value in redshift table")
        z_range = redshift[redshift <= zmax]
        comoving_volume = dVdz[redshift <= zmax] / (1. + z_range)

        if model == 'SH':
            z_pdf_comov = ff.SH(z_range, a=params['a'], b=params['b'], zm=params['zm'])
        elif model == 'BExp':
            z_pdf_comov = ff.BExp(z_range, a=params['a'], b=params['b'], zm=params['zm'])
        elif model == 'BPL':
            z_pdf_comov = ff.BPL_z(z_range, a=params['a'], b=params['b'], zm=params['zm'])
        elif model == 'D06':
            z_pdf_comov = ff.D06(z_range, a=params['a'], b=params['b'], c=params['c'], d=params['d'])
        elif model == 'qD06':
            z_pdf_comov = ff.qD06(z_range, SFR=params['SFR'])
        elif model == 'P16':
            z_pdf_comov = ff.P16(z_range, gamma_0=0.0204, gamma_1=1.8069, gamma_2=3.1724, gamma_3=7.2690)
        else:
            raise NotImplementedError

        pseudo_collapse_rate = integrate.trapz(z_pdf_comov/z_pdf_comov[0] * comoving_volume, z_range)
        self.normalization['pseudo_collapse_rate'] = pseudo_collapse_rate
        z_pdf_unnormed = z_pdf_comov * comoving_volume
        z_pdf = z_pdf_unnormed / z_pdf_unnormed.sum()

        z = np.random.choice(z_range, size=Nb_GRBs, p=z_pdf)

        self.properties['z'] = z
        self.properties['D_L'] = Lum_dist(z, cosmo)

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting z pdf")
            fig, axes = plt.subplots(2, figsize=(6, 7), tight_layout=True)
            delta_z = z_range[1]-z_range[0]
            axes[0].hist(z, bins=50, density=True, color='lightgray', label='Drawings',
                         edgecolor='k', linewidth=0.5)
            axes[0].plot(z_range, z_pdf/delta_z, label='PDF')
            axes[0].set_ylabel('z PDF')
            axes[0].set_yscale('log')
            axes[0].legend()
            axes[1].plot(z_range, z_pdf_comov, label='comoving rate')
            axes[1].set_ylabel(r'$\rm[yr^{-1}\,Mpc^{-3}]$')
            if model == 'P16':
                df = pd.read_csv('../catalogs/BAT6_cat/BAT6ext_GRB_formation_rate.txt',
                                 sep='\t', header=1, names=['1+z','GRB rate', 'err'])
                axes[1].errorbar(df['1+z']-1., df['GRB rate']/11, yerr=df['err']/11, fmt='.')
                axes[1].set_xlim(0,8)
                axes[1].set_ylabel('Arbitrary units')
            axes[1].legend()
            axes[1].set_xlabel('Redshift')
            axes[1].set_yscale('log')
            if savefig:
                self.save_fig(fig, 'z_draw.pdf')
        return z.copy()

    def draw_Ep(self, Nb_GRBs=None, model='LN', L=None, run_mode=None, savefig=False, **params):
        """
            Draw Ep distribution
        """

        self.parameters['peak_energy_distribution'] = {'model':model, **params}

        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        L0 = 1.6e52  # erg/s

        if model == 'Fixed':
            Ep = params["Ep0"]*np.ones(Nb_GRBs)

        elif model == 'LN':
            t = np.random.normal(np.log10(params["Ep0"]), params["sigmaEp"], Nb_GRBs)
            Ep = 10.**t

        elif model == 'A':
            t = np.random.normal(0, params["sigmaEp"], Nb_GRBs)
            if L is None:
                if 'L' not in self.properties.keys():
                    raise ValueError("Could not find L for Amati Ep drawing")
                L = self.properties['L']
            Ep = params["Ep0"] * (L/L0)**params["alpha_amati"] * 10.**(t * np.sqrt(1. + params["alpha_amati"]**2))

        else:
            raise ValueError("Invalid model for Ep drawing")

        self.properties['Ep'] = Ep

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting Ep pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            ax.hist(np.log10(Ep),
                    bins=np.arange(-1, 4, 0.1),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(Ep [keV])')
            ax.set_ylabel('Ep PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 'Ep_draw.pdf')
        return Ep.copy()

    def draw_spec(self, Nb_GRBs=None, model='FBand', run_mode=None, savefig=False, data_dir=root_dir/'data', **params):
        """
            Draw the spectral parameters alpha, beta and ktild
        """

        self.parameters['spectral_shape'] = {'model':model, **params}

        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if model == 'FBand':
            _ktild = f90f.f90f.calc_ktild(alpha=params['alpha'], beta=params['beta'])
            ktild = _ktild * np.ones(Nb_GRBs)
            alpha = params['alpha'] * np.ones(Nb_GRBs)
            beta = params['beta'] * np.ones(Nb_GRBs)

        elif (model == 'GBM_Band_old') or (model == 'GBM_Band'):
            if (model == 'GBM_Band_old'):
                alpha_file = 'alpha_GBM.txt'
                beta_file = 'beta_GBM.txt'
            elif (model == 'GBM_Band'):
                alpha_file = 'good_alpha_GBM.txt'
                beta_file = 'good_beta_GBM.txt'

            alpha = st.draw_from_cdf_file(data_dir/alpha_file, N_draws=Nb_GRBs)
            beta = st.draw_from_cdf_file(data_dir/beta_file, N_draws=Nb_GRBs)

            if (model == 'GBM_Band_old'):  # because of wrong convention
                alpha = -alpha
                beta = -beta
                beta = np.where((beta == 2), 2.01, beta)

            ktild = f90f.f90f.calc_ktild(alpha=alpha, beta=beta)

        else:
            raise ValueError("Invalid model for spec drawing")

        self.properties['alpha'] = alpha
        self.properties['beta'] = beta
        self.properties['ktild'] = ktild

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting spectral params pdf")
            fig, axar = plt.subplots(2, figsize=(6, 7), tight_layout=True)
            for i, item, name in zip([0, 1], [alpha, beta], ['alpha', 'beta']):
                axar[i].hist(item,
                             bins=30,
                             density=True, color='lightgray',
                             label='Drawings', edgecolor='k', linewidth=0.5)
                axar[i].set_xlabel(name)
                axar[i].set_ylabel(name + ' PDF')
                axar[i].legend()
            if savefig:
                self.save_fig(fig, 'alpha_beta_draw.pdf')
        return alpha.copy(), beta.copy(), ktild.copy()

    def draw_t90(self, z_med=None, Nb_GRBs=None, run_mode=None, savefig=False, **params):
        """
            Draw t90 from LogNormal distribution corrected by the median
            redshift of the population on which the t90obs distribution
            was adjusted.
        """

        self.parameters['t90obs_distribution'] = params

        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if z_med is None:
            if 'z_med' in params.keys():
                z_med = params['z_med']
            else:
                log.warning("In draw_t90, no z_med was found to correct the t90 distribution."
                            " A default value of z_med = 1.7 was used but you should check this.")
                z_med = 1.7

        # Correct the observed mu with the median redshift of the
        # sample on which it was derived
        mu_t90 = params['mu']-np.log10(1.+z_med)
        sigma_t90obs = params['sigma']
        t90 = 10.**(np.random.normal(mu_t90, sigma_t90obs, Nb_GRBs))

        self.properties['t90'] = t90

        if 'z' in self.properties.keys():
            self.properties['t90obs'] = t90 * (1. + self.properties['z'])
            t90obs_exists = True
        else:
            log.warning("No T90obs was calculated because 'z' could not be found")
            t90obs_exists = False

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting t90 pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            logt90min = mu_t90 - 5*sigma_t90obs
            logt90max = mu_t90 + 5*sigma_t90obs
            ax.hist(np.log10(t90),
                    bins=np.arange(logt90min, logt90max, 0.1),
                    density=True, color='lightgray', alpha=0.7,
                    label='T90 Drawings', edgecolor='k', linewidth=0.5)
            if t90obs_exists:
                z_med_pop = np.median(self.properties['z'])
                ax.hist(np.log10(self.properties['t90obs']),
                        bins=np.arange(logt90min+(np.log10(1+z_med_pop)),
                                       logt90max+(np.log10(1+z_med_pop)), 0.1),
                        density=True, color='darkgray', alpha=0.8, zorder=0,
                        label='T90obs Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(T90 [s])')
            ax.set_ylabel('T90 PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 't90_draw.pdf')
        return t90.copy()

    def draw_Cvar(self, Nb_GRBs=None, t90obs=None, run_mode=None, savefig=False, **params):
        """
            Draw Cvar from t90obs correlated distribution
        """

        self.parameters['Cvar_distribution'] = params

        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if t90obs is None:
            if 't90obs' not in self.properties.keys():
                raise ValueError("Could not find t90obs for Cvar drawing")
            t90obs = self.properties['t90obs']

        mu_Cvar = params['mu']
        sigma_Cvar = params['sigma']
        correl_slope = params['correl_slope']
        t = np.random.normal(mu_Cvar, sigma_Cvar, Nb_GRBs)
        Cvar = 10.**(t + correl_slope * np.log10(t90obs))
        Cvar = np.where((Cvar >= 1), np.ones(Cvar.shape), Cvar)

        self.properties['Cvar'] = Cvar

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting Cvar pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            ax.hist(np.log10(Cvar),
                    bins=np.arange(-2, 0.1, 0.05),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(Cvar)')
            ax.set_ylabel('Cvar PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 'Cvar_draw.pdf')
        return Cvar.copy()

    def draw_GRB_properties(self, cosmo, params, run_mode=None, savefig=False):
        """
            Draw the various GRB properties used in the GRB population
            code. Also calculate some others from the drawn properties.
        """

        self.draw_z(**params['redshift_distribution'], cosmo=cosmo, run_mode=run_mode, savefig=savefig)
        self.draw_L(**params['luminosity_function'], run_mode=run_mode, savefig=savefig)
        self.draw_Ep(**params['peak_energy_distribution'], run_mode=run_mode, savefig=savefig)
        self.draw_spec(**params['spectral_shape'], run_mode=run_mode, savefig=savefig)
        self.properties['Epobs'] = self.properties['Ep']/(1. + self.properties['z'])

        # The median redshift of the population of GBM_bright is needed
        # to calculate t90 (which is defined as pht_pflx in BATSE band
        # above 0.9 ph/s/cm2)
        self.calc_peak_photon_flux(instruments={'BATSE':{'Emin':50,'Emax':300}})
        GBM_bright_selection = self.properties['pht_pflx_BATSE'] >= 0.9
        z_med = np.median(self.properties[GBM_bright_selection]['z'])
        self.draw_t90(**params['t90obs_distribution'], z_med=z_med, run_mode=run_mode, savefig=savefig)
        self.draw_Cvar(**params['Cvar_distribution'], run_mode=run_mode, savefig=savefig)

        self.properties['Eiso'] = self.properties['L'] * self.properties['Cvar'] * self.properties['t90']

        return self.properties

    def draw_GRB_properties_for_MCMC(self, cosmo, params, instruments):
        """
            Draw only the core properties necessary to run the MCMC
            exploration. This method is optimized for speed so not all
            checks are performed to see if the right input is used.
        """

        self.draw_z(**params['redshift_distribution'], cosmo=cosmo)
        self.draw_L(**params['luminosity_function'])
        self.draw_Ep(**params['peak_energy_distribution'])
        self.draw_spec(**params['spectral_shape'])

        self.properties['Epobs'] = self.properties['Ep']/(1. + self.properties['z'])
        # Use the function from physics module to avoid unecessary
        # checks present in the class method
        ph.calc_peak_photon_flux(GRB_prop=self.properties,
                                 instruments=instruments)

        return self.properties

    def calculate_quantities(self, instruments, samples, ECLAIRs_prop=None):
        """
            Convenience function to calculate:
            - peak photon flux
            - peak energy flux
            - photon fluence
            - energy fluence
            And then the detection probability for the included samples
        """
        self.calc_peak_photon_flux(instruments, ECLAIRs_prop=ECLAIRs_prop)
        self.calc_peak_energy_flux(instruments, ECLAIRs_prop=ECLAIRs_prop)
        self.calc_photon_fluence(instruments)
        self.calc_energy_fluence(instruments)
        if ECLAIRs_prop is not None:
            self.calc_det_prob(samples, **ECLAIRs_prop)
        else:
            self.calc_det_prob(samples)
        return

    def calc_peak_photon_flux(self, instruments, ECLAIRs_prop=None):
        """
            Calculate the peak photon flux for every GRB in the
            population.
        """
        necessary_prop = ['L', 'z', 'Ep', 'alpha', 'beta', 'D_L']
        self._check_properties(necessary_prop, func_name='peak photon flux')
        # This is to avoid recalculation BATSE pflx if it's already
        # been calculated for the t90 drawings
        if 'pht_pflx_BATSE' in self.properties.columns:
            _instruments = instruments.copy()
            _instruments.pop('BATSE')
        else:
            _instruments = instruments
        ph.calc_peak_photon_flux(GRB_prop=self.properties,
                                 instruments=_instruments,
                                 ECLAIRs_prop=ECLAIRs_prop)
        return

    def calc_peak_energy_flux(self, instruments, ECLAIRs_prop=None):
        """
            Calculate the peak energy flux for every GRB in the
            population.
        """
        necessary_prop = ['L', 'z', 'Ep', 'alpha', 'beta', 'D_L']
        self._check_properties(necessary_prop, func_name='peak energy flux')
        ph.calc_peak_energy_flux(GRB_prop=self.properties,
                                 instruments=instruments,
                                 ECLAIRs_prop=ECLAIRs_prop)
        return

    def calc_photon_fluence(self, instruments):
        """
            Calculate the photon fluence in units of ph/cm2 over the T90
            of the burst.
        """
        necessary_prop = ['t90obs', 'Cvar']
        self._check_properties(necessary_prop, func_name='photon fluence')
        ph.calc_photon_fluence(GRB_prop=self.properties,
                               instruments=instruments)
        return

    def calc_energy_fluence(self, instruments):
        """
            Calculate the photon fluence in units of ph/cm2 over the T90
            of the burst.
        """
        necessary_prop = ['t90obs', 'Cvar']
        self._check_properties(necessary_prop, func_name='energy fluence')
        ph.calc_energy_fluence(GRB_prop=self.properties,
                               instruments=instruments)
        return

    def calc_det_prob(self, samples, **ECLAIRs_prop):
        """
            Calculates the detection probability for the included
            samples
        """
        ph.calc_det_prob(GRB_prop=self.properties,
                         samples=samples, **ECLAIRs_prop)
        return

    def _check_properties(self, necessary_prop, func_name):
        """
            Convenience, private function to check if the inputs
            required for a given calculation are present in the
            population properties
        """
        for prop in necessary_prop:
            if prop not in self.properties.keys():
                raise KeyError(f"Property {prop} must exist before you attempt to calculate the"
                               f" {func_name}.")

    def create_mock_constraints(self, constraints=None):
        """
            Create the mock constraints from the current population.
            obs_constraints is a dictionary with the necessary
            information about each constraint.
        """

        if constraints is None:
            with open(root_dir/'init/obs_constraints.yml', 'r') as f:
                constraints = yaml.safe_load(f)
            load_observational_constraints(constraints)

        for name, constraint in constraints.items():
            bins = constraint['bins']
            val_min = constraint['val_min']
            prop_min = constraint['prop_min']
            quantity = constraint['quantity']
            cond = self.properties[prop_min] >= val_min
            mod, _u = np.histogram(self.properties[cond][quantity], bins=bins)
            self.mock_constraints[name] = {'val_min':val_min,
                                           'prop_min':prop_min,
                                           'quantity':quantity,
                                           'bins':bins,
                                           'hist_unnormed':mod,
                                           'err_unnormed':np.sqrt(mod)}
        return

    def compare_to_observational_constraints(self, constraints, method='chi2'):
        """
            First normalize the population to the observational
            constraints then calculate the likelihood of the current
            population.
        """
        lnL_tot = 0.0
        chi2_tot = 0.0
        for name, constraint in constraints.items():
            # Normalize to observations
            model = self.mock_constraints[name]['hist_unnormed']
            error = self.mock_constraints[name]['err_unnormed']
            norm = self._normalize_to_constraint(mod=model,
                                                 obs=constraint['hist'],
                                                 err=constraint['err'])
            self.mock_constraints[name]['norm'] = norm
            self.mock_constraints[name]['hist'] = norm * model
            self.mock_constraints[name]['err'] = norm * error
            # Calculate likelihood
            if method == 'chi2':
                chi2 = st.chi2(mod=self.mock_constraints[name]['hist'],
                               obs=constraint['hist'],
                               err=constraint['err'])
                self.likelihood_params['_'.join(['chi2',name])] = chi2
                chi2_tot += chi2
                lnL = -0.5 * chi2
            elif method == 'pBIL':
                lnL = st.pBIL(mod=self.mock_constraints[name]['hist'],
                              obs=constraint['hist'],
                              sum_ln_oi_factorial=constraint['sum_ln_oi_factorial'])
                # If using pBIL, add a factor of 10 weight to the eBAT6
                # constraint to make it more impactful
                # if name == 'eBAT6':
                #     lnL *= 10
            self.likelihood_params['_'.join(['lnL',name])] = lnL
            lnL_tot += lnL

        self.likelihood_params['lnL_tot'] = lnL_tot
        if method == 'chi2':
            self.likelihood_params['chi2_tot'] = chi2_tot

        return

    def normalize_from(self, nGRB0, nGRB0_err=None):
        """
            Normalize the population to a given local LGRB comoving
            event rate nGRB0 [yr-1 Mpc-3]
        """
        R_intr = nGRB0 * self.normalization['pseudo_collapse_rate']
        T_sim = self.Nb_GRBs/R_intr
        if nGRB0_err is not None:
            R_intr_err = R_intr * nGRB0_err/nGRB0
            T_sim_err = T_sim * R_intr_err/R_intr
        else:
            R_intr_err = np.nan
            T_sim_err = np.nan

        self.normalization['nGRB0'] = nGRB0
        self.normalization['nGRB0_err'] = nGRB0_err
        self.normalization['R_intr'] = R_intr
        self.normalization['R_intr_err'] = R_intr_err
        self.normalization['T_sim'] = T_sim
        self.normalization['T_sim_err'] = T_sim_err
        return

    def normalize_to_Stern(self):
        """
            Normalize the population using the Stern constraints.
            This yields T_sim the duration of the simulation
        """
        if 'Stern' not in self.mock_constraints.keys():
            raise KeyError("You must first create_mock_constraints for Stern before trying to"
                           " normalize an LGRB population")
        N_BATSE = np.sum(self.mock_constraints['Stern']['hist_unnormed'])
        # Simulation duration
        T_sim = N_BATSE / R_tot_BATSE
        T_sim_err = T_sim * R_tot_BATSE_err/R_tot_BATSE
        # Intrinsic Rate
        R_intr = self.Nb_GRBs/T_sim  # LGRB/yr in 4 pi above Lmin
        R_intr_err = R_intr * T_sim_err/T_sim
        # Local LGRB comoving rate density
        nGRB0 = R_intr / self.normalization['pseudo_collapse_rate']
        nGRB0_err = nGRB0 * R_intr_err/R_intr

        self.normalization['T_sim'] = T_sim
        self.normalization['T_sim_err'] = T_sim_err
        self.normalization['R_intr'] = R_intr
        self.normalization['R_intr_err'] = R_intr_err
        self.normalization['nGRB0'] = nGRB0
        self.normalization['nGRB0_err'] = nGRB0_err

        # Second method: chi2 minimisation
        if 'norm' not in self.mock_constraints['Stern'].keys():
            log.warning('Could not calculate T_sim_from_chi2 since no normalization'
                        ' was found for Stern')
        else:
            norm = self.mock_constraints['Stern']['norm']
            T_sim_from_chi2 = T_live_BATSE/norm
            self.normalization['T_sim_from_chi2'] = T_sim_from_chi2

        return

    def _normalize_to_constraint(self, mod, obs, err):
        """
            Estimate the best normalization factor using chi2 minimization.
        """
        norm = np.sum(obs*mod/err**2) / np.sum((mod/err)**2)
        return norm

    def save_to(self, fname=None, save_memory=True):
        """
            Save the population to a pickle file.
            You can set save_memory to True to save memory storage; in
            this case the drawings for each LGRB and the associated
            quantities wont be saved.
        """
        if fname is None:
            fname = self.output_dir/'GRB_population'
        elif fname is not None and not isinstance(fname, Path):
            fname = Path(fname)

        if not fname.parent.exists():
            fname.parent.mkdir()

        if save_memory:
            gp_to_save = GRBPopulation(Nb_GRBs=0)
            for attr, val in self.__dict__.items():
                if attr != 'properties':
                    gp_to_save.__dict__[attr] = val
        else:
            gp_to_save = self

        with open(fname, 'wb') as f:
            pickle.dump(gp_to_save, f)
        log.info('Saved GRB population to {}'.format(fname))
        return fname

    def save_fig(self, fig, filename):
        try:
            fig.savefig(self.output_dir/filename)
        except FileNotFoundError:
            self.output_dir.mkdir()
            log.info(f"Created '{self.output_dir}' to save {filename} ")
            fig.savefig(self.output_dir/filename)
        return

    def summary(self, full_width=80, cell_width=None, sides='|', thick_line="=", thin_line='-'):
        """
            A short summary of the current LGRB population returned as
            a string.
        """
        if full_width % 2 != 0:
            full_width += 1
            log.warning("In GRBPopulation.summary(): the width you are asking for is odd."
                        " It was increased by 1 for esthetic purposes.")
        try:
            max_len = max(map(len, self.properties)) + 1
        except ValueError:
            max_len = 14

        if cell_width is None:
            cell_width = np.max([max_len, 14])
        center_width = full_width - 2*len(sides)
        half_width = int(center_width / 2)
        summary = "\n" + full_width * thick_line + "\n"
        summary += sides + "SUMMARY".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"

        # General configuration
        _Nb_GRBs_string = f"Nb_GRBs =".rjust(half_width)
        _Nb_GRBs = f" {self.Nb_GRBs:.2e}".ljust(half_width)
        _output_dir_string = f"Output directory =".rjust(half_width)
        _output_dir = f" {self.output_dir.stem}".ljust(half_width)
        summary += sides + _Nb_GRBs_string + _Nb_GRBs + sides + "\n"
        summary += sides + _output_dir_string + _output_dir + sides + "\n"
        summary += full_width * thick_line + "\n"

        # Parameters
        summary += sides + "Parameters".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        for distr_name, params in self.parameters.items():
            summary += sides + distr_name.center(center_width) + sides + "\n"
            summary += full_width * thin_line + "\n"
            for param, value in params.items():
                _key = f"{param} =".rjust(half_width)
                _val = f" {value}".ljust(half_width)
                _summary = sides + _key + _val + sides + "\n"
                summary += _summary
            summary += full_width * thin_line + "\n"

        # Properties
        summary += full_width * thick_line + "\n"
        summary += sides + "Properties".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        summary += '|'.join([' prop'.ljust(cell_width),
                             'median'.center(cell_width),
                             'stdev'.center(cell_width),
                             'min'.center(cell_width),
                             'max'.center(cell_width)]) + "\n"
        summary += full_width * thin_line + "\n"
        for key, prop in self.properties.items():
            med = np.median(prop)
            std = np.std(prop)
            _min = prop.min()
            _max = prop.max()
            _summary = '|'.join([f" {key}".ljust(cell_width),
                                 f"{med:< .4e}".center(cell_width),
                                 f"{std:< .4e}".center(cell_width),
                                 f"{_min:< .4e}".center(cell_width),
                                 f"{_max:< .4e}".center(cell_width)])
            summary += _summary.ljust(full_width) + "\n"

        # Likelihood
        summary += full_width * thick_line + "\n"
        summary += sides + "Likelihood".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        for key, val in sorted(self.likelihood_params.items()):
            _key = f"{key} =".rjust(half_width)
            _val = f" {val:.3f}".ljust(half_width)
            _summary = sides + _key + _val + sides + "\n"
            summary += _summary

        # Normalization
        summary += full_width * thick_line + "\n"
        summary += sides + "Normalization".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        for key, val in sorted(self.normalization.items()):
            _key = f"{key} =".rjust(half_width)
            _val = f" {val:.3e}".ljust(half_width)
            _summary = sides + _key + _val + sides + "\n"
            summary += _summary

        summary += full_width * thick_line + "\n"

        return summary
