/*
    dczheng
    created 2019-11-06
*/

#include "allvars.h"

hid_t h5_dsp, h5_attr, h5_ds;
void hdf5_write_attr_scalar(
        hid_t h5_obj, 
        hid_t h5_type,
        char *name,
        void *data
        )  {

    h5_dsp = H5Screate( H5S_SCALAR );
    h5_attr = H5Acreate( h5_obj,
                         name,
                         h5_type,
                         h5_dsp,
                         H5P_DEFAULT );
    H5Awrite( h5_attr, h5_type, data );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );
}

void hdf5_write_attr_nd( 
        hid_t h5_obj,
        hid_t h5_type,
        hsize_t *h5_dims,
        int ndims,
        char *name,
        void *data
        ) {

    h5_dsp = H5Screate_simple( ndims, h5_dims, NULL );
    h5_attr = H5Acreate( h5_obj,
                         name,
                         h5_type,
                         h5_dsp,
                         H5P_DEFAULT );
    H5Awrite( h5_attr, h5_type, data );
    H5Aclose( h5_attr );
    H5Sclose( h5_dsp );
}

void hdf5_write_data( 
        hid_t h5_obj,
        hid_t h5_type,
        hsize_t *h5_dims,
        int ndims,
        char *name,
        void *data
        ) {

    h5_dsp = H5Screate_simple( ndims, h5_dims, NULL );
    h5_ds = H5Dcreate( h5_obj,
                       name,
                       h5_type,
                       h5_dsp,
                       H5P_DEFAULT );
    H5Dwrite( h5_ds,
               h5_type,
               h5_dsp,
               H5S_ALL,
               H5P_DEFAULT,
               data );
    H5Sclose( h5_dsp );
    H5Dclose( h5_ds );

}
