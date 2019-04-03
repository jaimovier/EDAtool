# EDAtool
Matlab tool for processing Electrodermal Activity (EDA) signals

# MATLAB usage 
Output from 'help EDAtool' in command window

    [EDA_f, EDA_x, EDA_P, EDA_T, Q] = ...
      EDAtool(EDA, SR, debug, InterOption, ResampleOption, par)
 
  Returns low-pass version of the signal & removes artefacts, with the
  option of interpolating between the artefact. Phasic (fast changes) and 
  Tonic (slow changes) components are obtained via IIR filtering. The 
  function resamples the EDA signal to 50Hz, with the option of 
  resampling back to the original SR.
 
  EDA_f: low-pass filtered EDA
  EDA_x: EDA_f without artefacts (NaN or Interpolation)
  EDA_P: Phasic EDA or EDR (with NaN)
  EDA_T: Tonic EDA or EDL
  Q: Confidence of EDA signal (%)
 
  EDA: EDA signal
  SR: Signal's Sample Rate
  debug: 0 for no debugging, 1 to debug (includes plot)
  ResampleOption: 0 for no resampling, 1 to resample
  InterOption: Option to interpolate between artefacts (0:off, 1:on)
  par is a vector with [RangeMin RangeMax WindowSize Threshold NaNSize]
  being:
    RangeMin: Lowest value in EDA range (e.g. open circuit)
    RangeMax: Highest value in EDA range (e.g. closed circuit)
    WindowSize: Size of artefact window in seconds.
    Threshold: Threshold for detecting artefacts within window (0-100)
    NaNSize: Size of NaN replacement vector for each artefact.
    
    Default values:
        debug = 0
        ResampleOption = 0
        RangeMin = 150 (for BioControl sensor with 10bit ADC)
        RangeMax = 490 (for BioControl sensor)
        WindowSize = 0.2 [s]
        Threshold = 4%
        NaNSize = 1.5 [s]
 
  Note: Values outside or equal to range limits will be considered
  artefacts.

# Reference
If you are using this tool for academic purposes, please cite the following paper:

Jaimovich, Javier, and R. Benjamin Knapp. 2015. “Creating Biosignal Algorithms for Musical Applications from an Extensive Physiological Database.” In Proceedings of the 2015 Conference on New Interfaces for Musical Expression (NIME 2015). Baton Rouge, LA.
