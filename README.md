## Calculate Prochlorococcus production profiles using DWCPN

This is a wrapper around the DWCPN primary production model library (https://github.com/JamesDingle/dwcpn). Currently referencing the branch 'adding-prochloro' while I am finishing up the implementation.

## System Requirements

You will need the following installed:
- Rust (I recommend installing via [rustup](https://rustup.rs))
- NetCDF4 and HDF4 libraries
  - I have tested with:
    - NetCDF4 v4.8.1
    - HDF5 v1.12.1
  - For macOS these are both available via home brew `brew install netcdf` `brew install hdf5`
  - For linux these are available through the native package managers on the mainstream linux distributions

## Inputfile requirements
This wrapper requires input data to be in a specific format currently which is a netcdf4 file on a lat/lon grid. The requirements are as follows:
- dimensions: ['lat', 'lon', 'depth']
- input variables
  - `chlor_a` - chlorophyll-a concentration (in mg m^-3)
  - `par` - Photosynthetically available radiation
  - `bathymetry` - (in metres)
  - `mld` - mixed layer depth (metres)
  - `sigma` - Chlorophyll profile parameter
  - `zm` - Chlorophyll profile parameter
  - `rho` - Chlorophyll profile parameter
  - `PI_pmb` - PI light curve parameter
  - `PI_alpha` - PI light curve parameter
  - `Pro_Surf` - Surface cell count for prochlorococcus
  - `Pro_Max` - maximum cell count for prochlorococcus
- output variables - these will be overwritten with the output data
  - `pp` - [lat, lon] - daily integrated primary production
  - `i_star_mean` - [lat, lon] - (not sure, need to get a definition)
  - `euphotic_depth` - [lat, lon] - depth at which available light goes below 1% of the surface value
  - `pro_total` - [depth, lat, lon] - total prochlorococcus production
  - `pro_1` - [depth, lat, lon] - prochlorococcus production from 1st kind of prochlorococcus (jad: fix naming)
  - `pro_2` - [depth, lat, lon] - prochlorococcus production from 2nd kind of prochlorococcus (jad: fix naming)

    
## Compilation & Usage

1. Clone this repo to your local machine 
2. run `cargo build --release` just to check that it compiles
3. run it using either (replacing values in <> brackets):
   1. `cargo run --release -- --inputfile <inputfile.nc> --jday <julian day of year>`
   2. `./target/release/pp_lib_example --inputfile <inputfile.nc> --jday <julian day of year>`