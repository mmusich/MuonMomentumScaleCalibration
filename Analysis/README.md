# Layer-by-layer corrections

## Recipe
The scripts to be used are aggregategrads.py and globalfith5pynoreduction.py, the first takes as input files the nTuples produced by ResidualGlobalCorrectionMakerTwoTrackG4e.

The scripts currently use narf (this can be changed), so Wremnants is needed. You might also need to clone scipy-muCal

```
git clone https://github.com/bendavid/scipy-MuCal.git -b obsopt 
```

and then use the versions of the scripts above in erc-asymow/MuonMomentumScaleCalibration.
