"""
Read in the liquid water Lorenz-Mie look-up table and interpolate over an array of effective radius.

Uses scale invariance, with base values valid at 355 nm. So we need to scale to 1064 nm (CPL depol channel).

Publication: Chemyakin Eduard, Stamnes Snorre, Burton Sharon P., Liu Xu, Hostetler Chris, Ferrare Richard, Cairns Brian, and Dubovik Oleg (2021) Improved Lorenz-Mie Look-Up Table for Lidar and Polarimeter Retrievals. Frontiers Remote Sensing. 2:711106. doi: 10.3389/frsen.2021.711106

Data/Software URL: https://science-data.larc.nasa.gov/polarimetry/
"""


import os
import numpy as np
import xarray as xr
import h5py
from scipy.interpolate import RegularGridInterpolator
from datetime import datetime

def load_liquid_lut(mat_filename):
    print(f"Loading Liquid LUT: {mat_filename.split('/')[-1]}...")
    with h5py.File(mat_filename, 'r') as f:
        # Access the top group
        cloud_grp = f['CloudLUT']
        
        # 1. Read Dimensions
        # Note: Matlab HDF5 is column-major. We often need to transpose (.T) 
        # to get [Real, Imag, Radius] in Python (Row-major).
        real_cri = np.squeeze(cloud_grp['RealPartsCRI'][:])
        imag_cri = np.squeeze(cloud_grp['ImagPartsCRI'][:])
        radii = np.squeeze(cloud_grp['RadiiGridBins'][:])
        angles = np.squeeze(cloud_grp['ScatteringAngles'][:])
        
        # 2. Access DataRecord Group
        data_rec = cloud_grp['DataRecord']
        
        # 3. Read Specific Optical Variables
        # Cbsc = Backscattering Cross Section (Already includes P11(180) logic!)
        # Shape usually: [Radius, Imag, Real] in Matlab -> [Real, Imag, Radius] in Python
        # CHECK SHAPES: print(C_bsc.shape) to be sure.
        C_bsc = data_rec['Cbsc'][:]
        C_ext = data_rec['Cext'][:]
        print(f"C_bsc Shape: {C_bsc.shape}")
        print(f"C_ext Shape: {C_ext.shape}")
        print(
            f"Grid Sizes: Real={len(real_cri)}, Imag={len(imag_cri)}, Rad={len(radii)}, Ang={len(angles)}"
        )
        
        # If the shape is [dim1, dim2, dim3], and radii has length dim3,
        # then no transpose needed. If radii corresponds to dim1, transpose.
        # Standard guess: Transpose is usually needed for Matlab -> Numpy.
        C_bsc = C_bsc.T
        C_ext = C_ext.T
        
    return real_cri, imag_cri, radii, C_bsc, C_ext
    
def reconstruct_array(h5_file, dataset_name, expected_shape_2d):
    """
    Reconstructs dense numpy array from HDF5 object references.
    """
    print(f"Processing {dataset_name}...")
    
    # Get the array of references (1400, 1)
    ref_array = h5_file['CloudLUT']['DataRecord'][dataset_name] 
    
    # Check inner shape from the first reference
    first_ref = ref_array[0, 0]
    inner_data = h5_file[first_ref][:]
    
    # Handle (700, 1) vs (700,)
    inner_shape = inner_data.shape
    if len(inner_shape) == 2 and inner_shape[1] == 1:
        inner_shape = (inner_shape[0],)
        
    final_shape = expected_shape_2d + inner_shape
    
    # Pre-allocate output array
    out_arr = np.zeros(final_shape, dtype=np.float32)
    
    # Loop and fill
    total_refs = expected_shape_2d[0] * expected_shape_2d[1]
    
    for i in range(total_refs):
        # Map 1D index to 2D (Real, Imag) - Row Major
        r_idx = i // expected_shape_2d[1]
        i_idx = i % expected_shape_2d[1]
        
        ref = ref_array[i, 0]
        data = h5_file[ref][:]
        
        if data.shape != inner_shape:
            data = data.reshape(inner_shape)
            
        out_arr[r_idx, i_idx, ...] = data
        
    return out_arr

# --- CONFIGURATION ---
input_h5_path = (
    'CloudLUT_100_100_0.350_700bins_1400CRI_ln2rKr_Twomey.mat'
)

# Dimensions
n_real = 56
n_imag = 25
n_rad = 700
n_ang = 203

# Variables to process
vars_3d = ['Cext', 'Cbsc', 'Csca']       # Shape: (Real, Imag, Radius)
vars_4d = ['C11', 'C12', 'C33', 'C34']   # Shape: (Real, Imag, Radius, Angle)

# --- MAIN EXECUTION ---
with h5py.File(input_h5_path, 'r') as f_in:
    # 1. EXTRACT COORDINATES
    real_vals = np.squeeze(f_in['CloudLUT']['RealPartsCRI'][:])
    imag_vals = np.squeeze(f_in['CloudLUT']['ImagPartsCRI'][:])
    
    if 'RadiiGridBins' in f_in['CloudLUT']:
        rad_vals = np.squeeze(f_in['CloudLUT']['RadiiGridBins'][:])
    else:
        print("Warning: 'Radii' not found, using placeholder.")
        rad_vals = np.logspace(np.log10(1e-3), np.log10(500), n_rad)

    if 'ScatteringAngles' in f_in['CloudLUT']:
        ang_vals = np.squeeze(f_in['CloudLUT']['ScatteringAngles'][:])
    else:
        print("Warning: 'ScatteringAngles' not found, using placeholder.")
        ang_vals = np.linspace(0, 180, n_ang)

    # 2. BUILD DATA DICTIONARY
    data_vars = {}
    
    # Process 3D Variables
    for var_name in vars_3d:
        data_numpy = reconstruct_array(f_in, var_name, (n_real, n_imag))
        
        da = xr.DataArray(
            data_numpy,
            dims=['real', 'imag', 'radius'],
            coords={'real': real_vals, 'imag': imag_vals, 'radius': rad_vals},
            # FIX: Correct units to volume-normalized dimensions (um^-1)
            attrs={'units': 'um^-1', 'long_name': f'{var_name} Volume-Normalized Coefficient'}
        )
        data_vars[var_name] = da

    # Process 4D Variables
    for var_name in vars_4d:
        data_numpy = reconstruct_array(f_in, var_name, (n_real, n_imag))
        
        if data_numpy.ndim == 4:
            dims = ['real', 'imag', 'radius', 'angle']
            coords = {'real': real_vals, 'imag': imag_vals, 'radius': rad_vals, 'angle': ang_vals}
        else:
            dims = ['real', 'imag', 'radius']
            coords = {'real': real_vals, 'imag': imag_vals, 'radius': rad_vals}
            
        da = xr.DataArray(
            data_numpy,
            dims=dims,
            coords=coords,
            attrs={'units': 'um^-1 sr^-1', 'long_name': f'Scattering Matrix Element {var_name}'}
        )
        data_vars[var_name] = da

    # 3. CREATE DATASET
    ds = xr.Dataset(data_vars)
    
    # Add Global Attributes
    ds.attrs = {
        'description': 'Reconstructed Liquid Cloud LUT',
        'history': f'Created {datetime.now().isoformat()}',
        'source_file': input_h5_path,
        'wavelength_lut_um': 0.355
    }
    
# Save Dataset to file
output_nc_path = 'SIR_Liquid_LUT.nc'

# Define compression (zlib=True) for all variables
encoding = {var: {'zlib': True, 'complevel': 4} for var in ds.data_vars}

print(f"Saving to {output_nc_path}...")
ds.to_netcdf(output_nc_path, encoding=encoding)
print("Done!")