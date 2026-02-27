import yaml
import os
import json
import gzip
import pickle
import argparse
import dctools
import hist
import matplotlib.pyplot as plt
from dctools import plot_fit as plotter
from typing import Any, IO
import numpy as np
from tqdm import tqdm
class config_input:
    def __init__(self, cfg):
        self._cfg = cfg 
    
    def __getitem__(self, key):
        v = self._cfg[key]
        if isinstance(v, dict):
            return config_input(v)

    def __getattr__(self, k):
        try:
            v = self._cfg[k]
            if isinstance(v, dict):
                return config_input(v)
            return v
        except:
            return None
    def __iter__(self):
        return iter(self._cfg)



class config_loader(yaml.SafeLoader):
    """YAML Loader with `!include` constructor."""
    def __init__(self, stream: IO) -> None:
        """Initialise Loader."""
        try:
            self._root = os.path.split(stream.name)[0]
        except AttributeError:
            self._root = os.path.curdir
        super().__init__(stream)


def construct_include(loader: config_loader, node: yaml.Node) -> Any:
    """Include file referenced at node."""
    filename = os.path.abspath(os.path.join(loader._root, loader.construct_scalar(node)))
    extension = os.path.splitext(filename)[1].lstrip('.')

    with open(filename, 'r') as f:
        if extension in ('yaml', 'yml'):
            return yaml.load(f, config_loader)
        elif extension in ('json', ):
            return json.load(f)
        else:
            return ''.join(f.readlines())

yaml.add_constructor('!include', construct_include, config_loader)

def main():
    parser = argparse.ArgumentParser(description='The Creator of Combinators')
    parser.add_argument("-i"  , "--input"   , type=str , default="./config/input_UL_2018_timgad-vbs.yaml")
    parser.add_argument("-v"  , "--variable", type=str , default="nnscore")
    parser.add_argument("-y"  , "--era"     , type=str , default='2018')
    parser.add_argument("-c"  , "--channel" , nargs='+', type=str)
    parser.add_argument("-s"  , "--signal"  , nargs='+', type=str)
    parser.add_argument('-n'  , "--name"    , type=str , default='')
    parser.add_argument('-p'  , "--plot"    , action="store_true")
    parser.add_argument('--rebin', type=int, default=1, help='rebin')
    parser.add_argument("--bins", 
            type=lambda s: [float(item) for item in s.split(',')], 
            help='input a comma separated list. ex: --bins="-1.2,0,1.2"'
    )
    parser.add_argument('--blind', action='store_true', help='blinding the channel')
    parser.add_argument('--checksyst', action='store_true')
    parser.add_argument("-d" ,  '--dd', type=bool , default=False )

    options = parser.parse_args()
    config = dctools.read_config(options.input)

    print(f'making: {options.channel} : {options.variable} : {options.era}')

    if len(options.channel) == 1:
        options.channel = options.channel[0]
    
    # make datasets per prcess
    datasets = {}
    signal = ""

    if options.name=='':
        options.name == options.channel
        
        
    datasets:Dict = dict()
    for name in config.groups:
        histograms = dict(
            filter(
                lambda _n: _n[0] in config.groups[name].processes,
                config.boosthist.items()
            )
        )
        
        p = dctools.datagroup(
            histograms = histograms,
            ptype      = config.groups[name].type,
            observable = options.variable,
            name       = name,
            xsections  = config.xsections,
            channel    = options.channel,
            luminosity = config.luminosity.value,
            rebin      = options.rebin,
            era        = options.era
        )
        
        datasets[p.name] = p
        if p.ptype == "signal":
            signal = p.name


    if options.plot:
        _plot_channel = plotter.add_process_axis(datasets)
        pred = _plot_channel.project('process','systematic', options.variable)[:hist.loc('data'),:,:]
        data = _plot_channel[{'systematic':'nominal'}].project('process',options.variable)[hist.loc('data'),:] 

        plt.figure(figsize=(6,7))
        ax, bx = plotter.mcplot(
            pred[{'systematic':'nominal'}].stack('process'),
            data=None if options.blind else data, 
            syst=pred.stack('process'),
        )
        
        try:
            sig_ewk = _plot_channel[{'systematic':'nominal'}].project('process', variable)[hist.loc('VBSZZ2l2nu'),:]   
            sig_qcd = _plot_channel[{'systematic':'nominal'}].project('process', variable)[hist.loc('ZZ2l2nu'),:]   
            sig_ewk.plot(ax=ax, histtype='step', color='red')
            sig_qcd.plot(ax=ax, histtype='step', color='purple')
        except:
            pass
    
        ymax = np.max([line.get_ydata().max() for line in ax.lines if line.get_ydata().shape[0]>0])
        ymin = np.min([line.get_ydata().min() for line in ax.lines if line.get_ydata().shape[0]>0])
    
        ax.set_ylim(0.001, 100*ymax)
        ax.set_title(f"channel {options.channel}: {options.era}")

        ax.set_yscale('log')
        plt.savefig(f'plot-{options.channel}-{options.variable}-{options.era}.pdf')


    if options.checksyst:        
        _plot_channel = plotter.add_process_axis(datasets)
        pred = _plot_channel.project('process','systematic', options.variable)[:hist.loc('data'),:,:]
        data = _plot_channel[{'systematic':'nominal'}].project('process',options.variable)[hist.loc('data'),:] 
        plotter.check_systematic(
            pred[{'systematic':'nominal'}].stack('process'),
            syst=pred.stack('process'),
            plot_file_name=f'check-sys-{options.channel}-{options.era}'
        )
    with open('eft-names.dat') as eft_file:
        eftnames = [n.strip() for n in eft_file.readlines()]
    output_file = f"cards-eft_ZZ2l2nu/{options.era}_rate.txt"
    with open('eft-names.dat') as eft_file:
        eftnames = [n.strip() for n in eft_file.readlines()]
    def _scale_dd_uncertainty(shape_tuple, year_tag):
        """
        Rescale DD ratio up/down shapes so their integrals match the
        dedicated MET-bin uncertainties used when deriving the weights.
        """
        if not isinstance(shape_tuple, tuple) or len(shape_tuple) != 2:
            return shape_tuple
        up_hist, down_hist = shape_tuple
        if up_hist is None or down_hist is None:
            return shape_tuple

        year_str = str(year_tag)
        if '2016' in year_str:
            low_up, low_down = 1.086, 0.924
        elif year_str in ('2017', '2018'):
            low_up, low_down = 1.019, 0.983
        else:
            return shape_tuple

        def _apply(hist_obj, factor):
            return hist_obj * factor

        scaled_up = _apply(up_hist, low_up)
        scaled_down = _apply(down_hist, low_down)
        return (scaled_up, scaled_down)
    for eftn in tqdm(eftnames):
        card_name = options.channel+options.era+eftn

        card = dctools.datacard(
            name = signal if len(options.name)==0 else options.name,
            channel= card_name
        )
        card.shapes_headers()
        
        data_obs = datasets.get("data").get("nominal") #.view(flow=True)
        
        card.add_observation(data_obs)

        for _, p in datasets.items():
            # print(p.name,p.get("nominal").sum().value)
            if 'eft' in p.name:
                with gzip.open(f"/eos/user/h/hgao/ZZTo2L2Nu/PKL/aQGC-new/{options.era}/{eftn}.pkl.gz", 'rb') as f:
                    file_data = pickle.load(f)
      
                histograms = dict(
                    filter(
                        lambda _n: _n[0] in "ZZ2JTo2L2Nu2J_EWK_aQGC_TuneCP5_13TeV-madgraph-pythia8",
                        file_data.items()
                    )
                )

                p = dctools.datagroup(
                    histograms = histograms,
                    ptype      = config.groups["eft_ZZ2l2nu"].type,
                    observable = options.variable,
                    name       = "eft_ZZ2l2nu",
                    xsections  = config.xsections,
                    channel    = options.channel,
                    luminosity = config.luminosity.value,
                    rebin      = options.rebin,
                    era        = options.era
                )
                datasets[p.name] = p
                print(eftn,p.get("nominal").sum().value)
            if len(p.to_boost().shape) == 0 or p.get("nominal").sum().value == 0:
                print(f"--> histogram for the process {p.name} is empty !")
                continue
            if p.ptype=="data":
                continue
            if not card.add_nominal(p.name, p.get("nominal"), p.ptype): continue
            if 'eft' in p.name:
                with open(output_file, 'a') as fout:
                    fout.write(f"{eftn} {card.rates[7][1]}\n")
            year = options.era.replace('APV','')
            era_s = options.era.replace('APV','preVFP')
            
            if "DY" in p.name:
                dd_shape = p.get(f"dataDrivenDYRatio_{year}")
                # dd_shape = _scale_dd_uncertainty(dd_shape, year)
                card.add_shape_nuisance(p.name, f"CMS_SMP23001_DY_dd_uncert_{year}", dd_shape, symmetrise=False)
                # card.add_auto_stat()
                continue


            if 'WW' not in p.name and 'WZ' not in p.name and 'DY' not in p.name and 'Top' not in p.name:
                card.add_log_normal_lumi(p.name, f"lumi_{year}", config.luminosity.uncer)
                card.add_log_normal_lumi(p.name, f"lumi_13TeV_correlated", config.luminosity.uncer_correlated)
                if "16" not in year:
                    card.add_log_normal_lumi(p.name, f"lumi_13TeV_1718", config.luminosity.uncer_correlated1718)

            # interference between QCD and EWK
            card.add_log_normal(p.name, f"CMS_SMP23001_Interference_{options.era}", 1.0798)

            # HEM 15/16 
            if "18" in year:
                card.add_shape_nuisance(p.name, f"CMS_HEM_2018"  , p.get("HEM"), symmetrise=False)

            # scale factors / resolution
            card.add_shape_nuisance(p.name, f"CMS_res_e_{year}"  , p.get("ElectronEn"), symmetrise=True)
            card.add_shape_nuisance(p.name, f"CMS_scale_m"  , p.get("MuonRoc")   , symmetrise=True)
            
            card.add_shape_nuisance(p.name, f"CMS_SMP23001_lept_sf_{options.era}", p.get("LeptonSF")  , symmetrise=False)
            card.add_shape_nuisance(p.name, f"CMS_SMP23001_trig_sf_{options.era}", p.get("triggerSF") , symmetrise=False)

            # JES/JES and UEPS 
            # card.add_shape_nuisance(p.name, f"CMS_scale_j_{options.era}", p.get("JES"), symmetrise=False) 

            card.add_shape_nuisance(p.name, f"CMS_scale_j_Absolute_{year}"      , p.get(f"JES_Absolute{year}")      , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_BBEC1_{year}"         , p.get(f"JES_BBEC1{year}")         , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_EC2_{year}"           , p.get(f"JES_EC2{year}")           , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_HF_{year}"            , p.get(f"JES_HF{year}")            , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_RelativeSample_{year}", p.get(f"JES_RelativeSample{year}"), symmetrise=False) 

            card.add_shape_nuisance(p.name, f"CMS_scale_j_Absolute"   , p.get("JES_Absolute")   , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_BBEC1"      , p.get("JES_BBEC1")      , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_EC2"        , p.get("JES_EC2")        , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_HF"         , p.get("JES_HF")         , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_RelativeBal", p.get("JES_RelativeBal"), symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_scale_j_FlavorQCD"  , p.get("JES_FlavorQCD")  , symmetrise=False) 
            card.add_shape_nuisance(p.name, f"CMS_res_j_{era_s}", p.get("JER"), symmetrise=False)
            card.add_shape_nuisance(p.name, f"CMS_scale_met_unclustered_energy_{era_s}", p.get("UES"), symmetrise=False)
            
            # no correlated over era's
            card.add_shape_nuisance(p.name, f"ps_fsr", p.get("UEPS_FSR"), symmetrise=False)
            card.add_shape_nuisance(p.name, f"ps_isr", p.get("UEPS_ISR"), symmetrise=False)
            
            # b-tagging uncertainties
            # btag_sf_bc_2016APV, btag_sf_light_2016APV
            try:
                card.add_shape_nuisance(p.name, f"CMS_btag_fixedWP_incl_light_uncorrelated_{era_s}" , p.get(f"btag_sf_light_{options.era}"), symmetrise=True)
                card.add_shape_nuisance(p.name, f"CMS_btag_fixedWP_comb_bc_uncorrelated_{era_s}"  , p.get(f"btag_sf_bc_{options.era}")   , symmetrise=False)
            except:
                pass
            
            # # b-tagging uncertainties correlated over years
            card.add_shape_nuisance(p.name, "CMS_btag_fixedWP_comb_bc_correlated"  , p.get("btag_sf_bc_correlated")   , symmetrise=False)
            card.add_shape_nuisance(p.name, "CMS_btag_fixedWP_incl_light_correlated" , p.get("btag_sf_light_correlated"), symmetrise=True)
            card.add_shape_nuisance(p.name, f"CMS_l1_ecal_prefiring_{era_s}" , p.get("prefiring_weight"), symmetrise=True)
            # # other uncertainties
            card.add_shape_nuisance(p.name, f"CMS_pileup_{year}", p.get("pileup_weight"), symmetrise=False)

            #QCD scale, PDF and other theory uncertainty
            if 'gg' not in p.name:
                card.add_qcd_scales(
                        p.name, f"QCDscale_{p.name}", 
                        [p.get("QCDScale0w"), p.get("QCDScale1w"), p.get("QCDScale2w")]
            )
            
            # PDF uncertaintites / not working for the moment
            card.add_shape_nuisance(p.name, f"CMS_SMP23001_pdf_{p.name}"   , p.get("PDF_weight"), symmetrise=False)
            card.add_shape_nuisance(p.name, f"CMS_SMP23001_alphaS_{p.name}", p.get("aS_weight" ), symmetrise=False)        
            
            # Electroweak Corrections uncertainties
            if 'WZ' in p.name:
                card.add_shape_nuisance(p.name, "CMS_SMP23001_ewk_corr_WZ", p.get("kEW"), symmetrise=False)
            if ('ZZ' in p.name) and ('EWK' not in p.name):
                card.add_shape_nuisance(p.name, "CMS_SMP23001_ewk_corr_ZZ", p.get("kEW"), symmetrise=False)
                
            # define rates
            if p.name  in ["WW","Top"]:
                if "vbs-EM" in card_name:
                    card.add_rate_param(f"CMS_SMP23001_NormWW_{options.era}", "vbs-EM*", p.name)
                elif "SR" in card_name:
                    card.add_rate_param(f"CMS_SMP23001_NormWW_{options.era}", card_name+'*', p.name)
            
            # define rate 3L categoryel 
            elif p.name in ["WZ"]:
                if "vbs-3L" in card_name:
                    card.add_rate_param(f"CMS_SMP23001_NormWZ_{options.era}", "vbs-3L*", p.name)
                elif "SR" in card_name:
                    card.add_rate_param(f"CMS_SMP23001_NormWZ_{options.era}", card_name+'*', p.name)
            
            # define rate for DY category
            # elif p.name in ["DY"]:
            #     if "DY" in card_name:
            #         card.add_rate_param(f"NormDY_{options.era}", "vbs-DY*", p.name)
            #     elif "SR" in card_name:
            #         card.add_rate_param(f"NormDY_{options.era}", card_name+'*', p.name)
            
            card.add_auto_stat()

        # saving the datacard
        card.dump()

if __name__ == "__main__":
    main()

