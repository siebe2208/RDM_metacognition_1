"""
This file contains extra functions for ECG preprocessing pipelines 
use with ecg_preprocessing.py

Embodied Computation Group 

@author: Siebe Everaerts 
"""

from scipy.signal import find_peaks
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import neurokit2 as nk
# ============================================================================
# Trial Timing 
# ============================================================================
"""
These functions are used to detect trial timings for ECG experiments.

Parameters
----------
method: string 
    Either 'PD', 'stamps' or 'PD_check' does not support other ways.
    'PD_check' is a method that combines photodiode with stamps to validate the photodiode.

data: array-like 
    Raw photodiode signal for PD detection. Array with stamps for time stamps detection.

SR: float
    Sampling rate of the signal (Hz). Used to set the minimum distance between peaks in PD.

skip: list of int, optional 
    The indices of detected peaks to ignore/remove. Default is None.
    Only when method is 'PD' or 'PD_check'.

plot: bool, optional
    Set to True to return a plot. Only when method is 'PD' or 'PD_check'.

data_stamps: array-like 
    Extra array with time stamps for 'PD_check'. Should be provided only in 'PD_check'.

trans_samples: bool
    If data_stamps is in sec rather than samples, set this to True. 
    Should be provided only in 'PD_check' when data_stamps is in seconds.

Returns
result : varies
    - If method is 'PD': returns `timing`, a NumPy array of detected photodiode peaks.
    - If method is 'stamps': returns `sample_ID`, a NumPy array of sample indices converted from timestamps.
    - If method is 'PD_check': returns a dictionary with:
        {
            'timing': timing,  # NumPy array of detected PD peaks
            'check': check     # correlation coefficient between PD peaks and stamps
        }
----------


"""

def trial_timing(method, data, SR , data_stamps = None, skip = None, plot = False, trans_samples = False):

    # Check if method is actually supported and if so, which one is provided
    valid = ["PD", "stamps", "PD_check"]
    if method not in valid:
        raise ValueError("the trial timing detection only supports photodiode or time stamps. "
                        "Did you mean 'PD' or 'stamps'?")

    if method == "PD": 
        result = PD_detection(data, SR, skip, plot)

    elif method == "stamps":
         result = stamp_detection(data, SR)

    elif method == "PD_check":
        timing = PD_detection(data, SR, skip, plot)
        check = PD_stamp_check(timing, data_stamps, SR, trans_samples, plot)
        result = {"timing": timing, "check": check}
        
    
    return result
        

def PD_detection(PD_data, SR, skip = None, plot = False):
   
    # Ensure the data is a NumPy array

    PD_data = np.asarray(PD_data)

    # Normalize the signal (z-score) for consistent peak detection across recordings

    Z_PD = (PD_data - np.mean(PD_data)) / np.std(PD_data)

    # Set peak detection threshold relative to the max of the normalized signal

    thresh = max(Z_PD) - 1

    # Detect peaks
    # - distance ensures peaks are at least SR/10 samples apart
    # - height sets the minimum peak amplitude

    PD_peaks, _ = find_peaks(Z_PD, distance=SR/10, height=thresh)

    # Optional: remove peaks specified in `skip`

    if skip is not None:
        mask = np.ones(len(PD_peaks), dtype=bool)
        mask[skip] = False
        PD_peaks = PD_peaks[mask]
    
    # Optional: plot the outcome
    if plot:
        nk.events_plot(PD_peaks, PD_data)
        plt.show()

    return PD_peaks


def stamp_detection(stamps, SR):

    # Ensure the input is a NumPy array
    stamps = np.asarray(stamps)
    
    # Convert timestamps (seconds) to sample indices
    sample_ID = stamps * SR

    return sample_ID


def PD_stamp_check(PD_peaks, stamps, SR = None, trans_samples = False, plot = False):

    # Compute the difference between consecutive PD peaks and stamp events
    PD_peak_diff = np.diff(PD_peaks)
    stamp_diff = np.diff(stamps)

    # If the timestamps are in seconds, convert to samples
    if trans_samples:
        stamp_diff = stamp_diff * SR

    # Calculate correlation between PD peak intervals and stamp intervals
    try:
        corr_matrix = np.corrcoef(PD_peak_diff, stamp_diff)
        corr_coef = corr_matrix[0, 1]
    except Exception as e:
        print("An error occurred while computing correlation:", e)
        corr_coef = None
    finally:
        # Print array lengths for debugging
        print("PD array length is:", len(PD_peak_diff))
        print("stamp array length is:", len(stamp_diff))

    # Optional: plot the differences for visual inspection
    if plot:
        try:
            plt.scatter(PD_peak_diff, stamp_diff, color="blue", marker="o")
            plt.show()  
        except Exception:
            # Silently ignore plotting errors, same errors already dealt with when calculating corr_coef
            pass

    return corr_coef
        






#data_behavior = pd.read_csv(r"C:/Users/User/OneDrive/Documenten/psychologie 1e Master/internship/data behavior/0263HRD/0263HRD_final.txt")
#data_behavior_BH = pd.read_csv(r"C:/Users/User/OneDrive/Documenten/psychologie 1e Master/internship/data behavior/0263HRD/0263HRD_BH_final.txt")
#data_combined = pd.concat([data_behavior, data_behavior_BH], ignore_index = True)
#trial_start = data_combined["StartListening"]


#data = pd.read_csv( r"C:\Users\User\OneDrive\Documenten\psychologie 1e Master\internship\unprocessed data\sub_0263_HRD_plux_2021-05-31_09_50.csv")
#PD_data = data["Photo_diode"].to_numpy()

#y = trial_timing("PD_check", PD_data, 1000, data_stamps = trial_start, plot = True, skip = list(range(120, 128)) )
#print(y["timing"])
#print(y["check"])
