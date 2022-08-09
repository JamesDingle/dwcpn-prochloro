// use dwcpn::{self, dwcpn::{dwcpn::{ModelSettings, ModelInputs, calc_pp}, modules::pp_profile::{calculate_ay, calculate_bw, calculate_bbr}}};


use dwcpn::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};
use dwcpn::{ModelInputs, ModelSettings, ProchloroInputs};
use dwcpn::dwcpn::dwcpn::calc_production;
use dwcpn::dwcpn::modules::config::DEPTH_PROFILE_COUNT;
use indicatif::{ProgressBar, ProgressStyle};
use clap::{Command, Arg};

use ndarray::{Array, Array3, s};
use netcdf;
use netcdf::Variable;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    // Argument parsing                                                                                                                                                                                                                         
    let args = Command::new("DWCPN Primary Production Model")
        .version("0.1.0")
        .arg(Arg::new("inputfile")
            .short('i')
            .long("inputfile")
            .help("location of netcdf file to run the model on")
            .required(true)
            .takes_value(true)
        )
        .arg(Arg::new("jday")
            .short('j')
            .long("jday")
            .help("Julian day of year (a.k.a. Ordinal)")
            .required(true)
            .takes_value(true)
        )
        .arg(Arg::new("mld")
            .short('m')
            .long("mixed-layer-depth")
            .help("If enabled, limits depth profile to mixed layer depth with surface chlor_a concentration propagated throughout")
            .required(false)
            .takes_value(false)
        )
        .arg(Arg::new("iom")
            .short('p')
            .long("iom-only")
            .help("If enabled, limit the computation to the moment iom (noon par maximum) is calculated")
            .required(false)
            .takes_value(false)
        )
        .get_matches();

    let filename = args.value_of("inputfile").unwrap();
    let jday = args.value_of_t::<u16>("jday").unwrap();


    // pre calculate bw/bbr/ay arrays for use in the pp_profile calculation later
    // this only needs to be done once for all pixels so I have it here, before we start
    // looping over pixels
    let bw = calculate_bw();
    let bbr = calculate_bbr();
    let ay = calculate_ay();

    println!("Processing file: {}", filename);

    let mut ncfile = netcdf::append(&filename)?;

    // TODO: write a nice HashMap function to return all of these from a list of preset variable names
    // instead of two lines for each
    let lat = &ncfile.variable("lat").unwrap();
    let lon = &ncfile.variable("lon").unwrap();
    let chl = &ncfile.variable("chlor_a").unwrap();
    let par = &ncfile.variable("par").unwrap();
    let bathymetry = &ncfile.variable("bathymetry").unwrap();
    let pi_alpha = &ncfile.variable("PI_alpha").unwrap();
    let pi_pmb = &ncfile.variable("PI_pmb").unwrap();
    let zm = &ncfile.variable("zm").unwrap();
    let mld = &ncfile.variable("mld").unwrap();
    let rho = &ncfile.variable("rho").unwrap();
    let sigma = &ncfile.variable("sigma").unwrap();
    let pro_surf = &ncfile.variable("Pro_Surf").unwrap();
    let pro_max = &ncfile.variable("Pro_Max").unwrap();

    let lat_data = lat.values::<f64>(None, None)?;
    let lon_data = lon.values::<f64>(None, None)?;
    let chl_data = chl.values::<f64>(None, None)?;
    let par_data = par.values::<f64>(None, None)?;
    let bathymetry_data = bathymetry.values::<f64>(None, None)?;
    let pi_alpha_data = pi_alpha.values::<f64>(None, None)?;
    let pi_pmb_data = pi_pmb.values::<f64>(None, None)?;
    let zm_data = zm.values::<f64>(None, None)?;
    let mld_data = mld.values::<f64>(None, None)?;
    let rho_data = rho.values::<f64>(None, None)?;
    let sigma_data = sigma.values::<f64>(None, None)?;
    let pro_surf_data = pro_surf.values::<f64>(None, None)?;
    let pro_max_data = pro_max.values::<f64>(None, None)?;


    const F64_FILLVALUE: f64 = 9969209968386869000000000000000000000.0;
    let mut pp_data = vec![F64_FILLVALUE; lat.len() * lon.len()];
    let mut spectral_i_star_mean_data = vec![F64_FILLVALUE; lat.len() * lon.len()];
    let mut euphotic_depth_data = vec![F64_FILLVALUE; lat.len() * lon.len()];

    // 3d prochloro output depth profiles
    let pro_profile_shape = (DEPTH_PROFILE_COUNT, lat.len(), lon.len());
    let mut pp_prochloro_profile: Array3<f64> = Array3::from_elem(pro_profile_shape, F64_FILLVALUE);
    let mut pro_total_profile: Array3<f64> = Array3::from_elem(pro_profile_shape, F64_FILLVALUE);
    let mut pro_1_profile: Array3<f64> = Array3::from_elem(pro_profile_shape, F64_FILLVALUE);
    let mut pro_2_profile: Array3<f64> = Array3::from_elem(pro_profile_shape, F64_FILLVALUE);

    let count = lat.len() * lon.len();
    let pb = ProgressBar::new(count as u64);
    pb.set_draw_delta(1000);
    pb.set_style(ProgressStyle::default_bar()
        .template("{msg} {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {percent}% [{pos:>7}/{len:7} @ {per_sec}] (ETA: {eta})")
        .progress_chars("#>-"));

    for y in 0..lat.len() {
        for x in 0..lon.len() {
            pb.inc(1);
            if check_if_filled(chl_data[[y,x]], chl)
                || check_if_filled(par_data[[y, x]], par)
                || check_if_filled(bathymetry_data[[y, x]], bathymetry)
                || check_if_filled(mld_data[[y, x]], mld)
                || check_if_filled(pi_alpha_data[[y, x]], pi_alpha)
                || check_if_filled(pi_pmb_data[[y, x]], pi_pmb)
                || check_if_filled(zm_data[[y, x]], zm)
                || check_if_filled(rho_data[[y, x]], rho)
                || check_if_filled(sigma_data[[y, x]], sigma)
                || check_if_filled(pro_max_data[[y, x]], sigma)
                || check_if_filled(pro_surf_data[[y, x]], sigma)
            {
                continue;
            }

            let inputs = ModelInputs {
                lat: lat_data[y],
                lon: lon_data[x],
                z_bottom: bathymetry_data[[y, x]],
                iday: jday,
                // iday: 165,
                alpha_b: pi_alpha_data[[y, x]],
                pmb: pi_pmb_data[[y, x]],
                z_m: zm_data[[y, x]],
                mld: mld_data[[y, x]],
                chl: chl_data[[y, x]],
                rho: rho_data[[y, x]],
                sigma: sigma_data[[y, x]],
                cloud: 0.0,
                yel_sub: 500.0,
                par: par_data[[y, x]],
                bw,
                bbr,
                ay,
            };

            let settings = ModelSettings {
                mld_only: args.is_present("mld"),
                iom_only: args.is_present("iom"),
                prochloro_inputs: Some(ProchloroInputs{
                    prochloro_surface: pro_surf_data[[y,x]],
                    prochloro_maximum: pro_max_data[[y,x]]
                })
            };

            match calc_production(&inputs, &settings) {
                Ok(model_output) => {

                    // write out single value results
                    pp_data[y * lon.len() + x] = model_output.pp_day.unwrap();
                    euphotic_depth_data[y * lon.len() + x] = model_output.euphotic_depth.unwrap();
                    spectral_i_star_mean_data[y * lon.len() + x] = model_output.par_noon_max.unwrap();

                    // write out the depth profile results
                    let pro_total_array = Array::from_vec(
                        model_output.pro_total_profile.unwrap().clone().to_vec()
                    );
                    pro_total_profile.slice_mut(s![.., y, x]).assign(&pro_total_array);

                    let pro_1_array = Array::from_vec(
                        model_output.pro_1_profile.unwrap().clone().to_vec()
                    );
                    pro_1_profile.slice_mut(s![.., y, x]).assign(&pro_1_array);

                    let pro_2_array = Array::from_vec(
                        model_output.pro_2_profile.unwrap().clone().to_vec()
                    );
                    pro_2_profile.slice_mut(s![.., y, x]).assign(&pro_2_array);

                    let pp_prochloro_array = Array::from_vec(
                        model_output.pp_prochloro_profile.unwrap().clone().to_vec()
                    );
                    pp_prochloro_profile.slice_mut(s![.., y, x]).assign(&pp_prochloro_array);

                },
                Err(e) => {
                    println!("{:?}", e);
                }
            }
        } // x loop
    } // y loop

    pb.finish_with_message("Complete!");

    let mut pp_var = ncfile.variable_mut("pp").unwrap();
    pp_var.put_values(&pp_data, None, None).expect("Could not write out PP values");

    let mut euph_var = ncfile.variable_mut("euphotic_depth").unwrap();
    euph_var.put_values(&euphotic_depth_data, None, None).expect("Could not write out Euphotic Depth values");

    let mut spectral_i_star_var = ncfile.variable_mut("i_star_mean").unwrap();
    spectral_i_star_var.put_values(&spectral_i_star_mean_data, None, None).expect("Could not write out i_star_mean values");

    let mut pro_total_var = ncfile.variable_mut("pro_total").unwrap();
    let pro_total_output = pro_total_profile.to_owned().into_raw_vec();
    pro_total_var.put_values(&pro_total_output, None, None).expect("Could not write out Pro Total values");

    let mut pro_1_var = ncfile.variable_mut("pro_1").unwrap();
    let pro_1_output = pro_1_profile.to_owned().into_raw_vec();
    pro_1_var.put_values(&pro_1_output, None, None).expect("Could not write out Pro Total values");

    let mut pro_2_var = ncfile.variable_mut("pro_2").unwrap();
    let pro_2_output = pro_2_profile.to_owned().into_raw_vec();
    pro_2_var.put_values(&pro_2_output, None, None).expect("Could not write out Pro Total values");

    let mut pp_prochloro_var = ncfile.variable_mut("prochlorococcus").unwrap();
    let pp_prochloro_output = pp_prochloro_profile.to_owned().into_raw_vec();
    pp_prochloro_var.put_values(&pp_prochloro_output, None, None).expect("Could not write out Pro Total values");

    Ok(())
}
fn check_if_filled<T: PartialEq + netcdf::Numeric>(value: T, variable: &Variable) -> bool {
    match variable.fill_value::<T>() {
        Ok(f) => match f {
            Some(fill_value) => fill_value == value,
            None => false
        },
        Err(_) => false
    }
}
