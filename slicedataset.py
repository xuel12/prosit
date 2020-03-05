# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def slice_hdf5(path, n_samples=None):
    import h5py
    import numpy as np
    import random  

    # Get a list of the keys for the datasets
    random.seed(100)
    with h5py.File(path, 'r') as f:
        with h5py.File('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/data/traintest_hcd_10k.hdf5', 'w') as f2:
            print(f.keys())
            dataset_list = list(f.keys())
            
#            fulllen = f[dataset_list[0]].shape[0]
#            idx = sorted(random.sample(range(0,fulllen), int(fulllen*prop)))
            idx = range(n_samples)
            
            for dset_name in dataset_list:
                print(dset_name)
#                print(f[dset_name].dtype)
#                print(len(f[dset_name].shape))
#                print(f[dset_name].shape)
                
                if (len(f[dset_name].shape) == 2):
#                    print(f[dset_name][:1, :])
                    dset = f2.create_dataset(dset_name, shape=(len(idx), f[dset_name].shape[1]), dtype=f[dset_name].dtype)
                    dset.write_direct(f[dset_name][idx, :])
                elif (len(f[dset_name].shape) == 3):
#                    print(f[dset_name][:1, :, :])
                    dset = f2.create_dataset(dset_name, shape=(len(idx), f[dset_name].shape[1], f[dset_name].shape[2]), dtype=f[dset_name].dtype)
                    dset.write_direct(f[dset_name][idx, :, :])
                else:
#                    print(f[dset_name][:1, ])
                    dset = f2.create_dataset(dset_name, shape=(len(idx), ), dtype=f[dset_name].dtype)
                    dset.write_direct(f[dset_name][idx, ])
            f2.close()
        f.close()
            
#        sequence_integer_dset = f['sequence_integer']
#        print(sequence_integer_dset.shape)
#        print(sequence_integer_dset.dtype)
#
#        ce_dset = f['collision_energy_aligned_normed']
#        print(ce_dset.shape)
#        charge_dset = f['precursor_charge_onehot']
#        print(charge_dset.shape)
    return dataset_list

slice_hdf5('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/data/traintest_hcd.hdf5', n_samples = 10000)


def read_hdf5(path, n_samples=None):
    import h5py
    
    # Get a list of the keys for the datasets
    with h5py.File(path, 'r') as f:
        print(f.keys())
        dataset_list = list(f.keys())
        for dset_name in dataset_list:
            print(dset_name)
            print(f[dset_name].dtype)
            print(f[dset_name].shape)
        f.close()

    return dataset_list

read_hdf5('/Users/xuel12/Documents/MSdatascience/CS7180AI/project/data/traintest_hcd_10k.hdf5')

