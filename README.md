# LK_toolbox (a.k.a gipseeg)
Personal Matlab Toolbox written for EEG processing and classification during the begginning of my thesis @GIPSA-lab (2013-2015). Not well maintained but released for educational purpose.

It is mostly unfinished job. Modules of my work has been released into other repositories in more advanced state. After 2015, the code was maintained in a different (private) repo.

If you found this repo and plan to use it (or part of it), please read below before.

### few words

Because I came from signal processing, I wasn't fond of the complex eeglab/fieldtrip structures with never ending attributes and parameters. So I tried to build the "minimum-viable" structure to process and classify EEG on my own and it was probably not better than EEG/fieldtrip structure anyway. While I tried to maintain this code at some point, the effort were dropped because only few engineers and PhD students were using it and I decided to share only piece of cleaner code in different repo. I finally, after many years of this code roting on my computer, decided to upload the ugly monster. Let say it is a monument of "how to not start a thesis".


### Install
to install LK TOOLBOX and related functions please, launch the script "installer" (run in matlab)
thus you'll have acces to all the functions and script. For some functions, external requirements be needed (see external section below)

# The toolbox

## The content
Of course the library is far from being well organized because of poor decisions during my thesis so don't expect to find everything perfectly in place. Also, because I had to work on several papers at the same time, sometimes functions have several versions co-existing to not loose compatibility of my code.

### LIB folder
This folder contains _mostly_ original work.

- preproc: preprocessing data, epoching, convertion of trigger channel, some basic math fonctions
- Epoch: averaging, and other methods for organizing and processing EEG epochs
- processing: preprocess your eeg, withen your data, cospectra, etc.
- fileio (save-load-convert): All operations of loading. Out toolbox should be fully compatible for loading, converting on the following format : .mat
- classif: classification methods using Riemannian Geometry
- hyperscanning: classification and processing of simultaneous recording of multiple brains
- metrics: more Riemannian Geometry metrics used that were not in other toolbox
- synthetic: code to simulate EEG
- CSTP: dev of the repo https://github.com/lkorczowski/acstp (please use that one instead)
- BSS: dev prior of the release of https://github.com/lkorczowski/CAJD (please use that one instead)
and many more...

### script 

A messy folder containing plenty of experiments. Not very worthy if you don't have the data in hand but usefull if you want to see how the functions in lib folder were used.

### external

I give some of the code that is mandatory for running my toolbox. But I'll probably need to add the two following folders if some functions doesn't work:
- eeglab https://sccn.ucsd.edu/eeglab/index.php
- fieldtrip http://www.fieldtriptoolbox.org/


## the gipseeg structure
While I tried that most of the funtion can work directly with matrix (and not complicated structures), it may be easier to use to encapsulate your data into a structure, e.g. for using `preprocessingEEG` or `plotEEG`.


### Dimensions and nomenclature :
- variable - name - (scale) - Fieldtrip : dimord (variable name)
- t - samples (T) FT: time (time)
- n - channels (N) FT: chan (label)
- k - trials (K) FT: rpt (trial)
- m - subjects (M) FT: subj (‘’)
- f - frequencies (F) FT: freq (freq)
- z - classes (Z) FT: rpt (trialinfo)
- c - conditions (C) FT: rpt (trialinfo)
_constrains T>N, T>M, K>>M_


| FIELDNAME       | SIZE/TYPE                   | DESCRIPTION                                                                           | fieldtrip  | eeglab  | autre |
|-----------------|-----------------------------|---------------------------------------------------------------------------------------|------------|---------|-------|
| Fs              | scalar                      | sample rate in Hz                                                                     | fsample    | ‘srate’ | fs    |
| Trigger         | [nb samples x1 ]            | Trigger channel of '0' with '1' at the start of each sweep. There are [nb epochs] '1' | sampleinfo |         |       |
| EpochClass      | [nb epochs x1]              | class of the sweeps (0 for Non-TARGET, 1 for TARGET).                                 | trialinfo  |         |       |
| Channels        | [nb samples x nb channels]  | continuous eeg recording                                                              |            |         |       |
| NoiseTrigger*   | [nb samples x1 ]            | the equivalent of Trigger but for the sweep of the noise.                             |            |         |       |
|                 |                             |                   By default, it takes the same.                                      |            |         |       |
| ElectrodesName* | {1 x nb channels}           | the names of the electrodes                                                           | label      |         |       |


ACSTPoptions is a structure with
```
     Epoch_size: scalar, the length of the epoch window (in samples)
LatencyCorr_max: scalar, the maximum of samples for the latency
                  correction. Set 0 to disable the Latency correction.
Mask_Electrodes: vector of the selectionned electrodes. Usefull for the
                 automatic subspace selection (BestPz) and latency
                 correction (Latency).
      Mask_Time: vector of the selectionned sample. Usefull for the
                 automatic subspace selection (BestPz) and latency
                 correction (Latency).
MaxIterLatency*: scalar, the maximum iteration allowed to compute the
                  latency correction. Default: 10.
   SubspaceDim*: vector, containing all the subspace dimension (nb electrodes)
                 to test in descent order.
                  By default, it is equal to (nb_channels:-1:(nb_channels/2))
computeClassLat*: vector, containing all the class tag in which you want
                  to compute the latency correction. By default, it
                  computes it for all classes but it could be long (for
                  instance you can skip the non-target).
       Weights*: Default: true.
                 option1(given) [nb epochs x1] vector, containing the weights for each
                  epoch if it is computed from an external function.
                 option2 (true/false) boolean. If true (default) the ACSTP compute the
                  weight for each epoch. If false, the weights are
                  desactivated (i.e. set to 1).
       DISPLAY*: Boolean, should the comparative result between the arithmetic ensemble
                  average and the ACSTP should be display at the end. Default: true.
```

ACSTPstruct is a structure with
```
             EA: the ensemble average before ACSTP
         EAcstp: the ensemble average corrected with latencies, weighted
                  and filtered + with the effect of overlapping
                  correction
    As Bs Bt At: such as Xhat(:,:,k)=As{indClass}*Bs{indClass}'*W(k)*X(:,:,k)*Bt{indClass}*At{indClass}'
                  Each filter is a cell containing a matrix filter for
                  each class
          Class: The tag and the order in which the filter are sorted
         BestPz: Orders of the best subspace for the ACSTP for the given
                 Class
        Weights: [nb epochs x1] the weights of each epoch
        Latency: [nb epochs x1] the corrected offset of each epoch
     Epoch_size: scalar, the length of the epoch window (in samples)

```

# associated publications
- L. Korczowski, M. Congedo and C. Jutten, "Single-trial classification of multi-user P300-based Brain-Computer Interface using riemannian geometry," 2015 37th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC), Milan, 2015, pp. 1769-1772.
-  Louis Korczowski. Méthodes pour l'électroencéphalographie multi-sujet et application aux interfaces cerveau-ordinateur. Traitement du signal et de l'image [eess.SP]. Université Grenoble Alpes, 2018. Français. ⟨NNT : 2018GREAT078⟩. ⟨tel-01961434v4⟩
-  Marco Congedo, Louis Korczowski, Arnaud Delorme, Fernando Lopes da Silva. Spatio-temporal common pattern: A companion method for ERP analysis in the time domain. Journal of Neuroscience Methods, Elsevier, 2016, 267, pp.74-88. ⟨10.1016/j.jneumeth.2016.04.008⟩. ⟨hal-01343026⟩
-  Louis Korczowski, Florent Bouchard, Christian Jutten, Marco Congedo. Mining the Bilinear Structure of Data with Approximate Joint Diagonalization. 24th European Signal Processing Conference (EUSIPCO 2016), EURASIP, Aug 2016, Budapest, Hungary. pp.667-671. ⟨hal-01357245⟩
-  Florent Bouchard, Louis Korczowski, Jérôme Malick, Marco Congedo. Approximate Joint Diagonalization within the Riemannian Geometry Framework. 24th European Signal Processing Conference (EUSIPCO 2016), Aug 2016, Budapest, Hungary. pp.210-214, ⟨10.1109/EUSIPCO.2016.7760240⟩. ⟨hal-01370052
- Louis Korczowski, Alexandre Barachant, Anton Andreev, Christian Jutten, Marco Congedo. ”BrainInvaders 2” :  an open source Plug & Play multi-user BCI videogame.   6th International Brain-Computer Interface Meeting (BCI Meeting 2016), BCI Society, May 2016, Pacific Grove, CA, UnitedStates. pp.10.3217/978-3-85125-467-9, ￿10.3217/978-3-85125-467-9-224￿. ￿hal-
