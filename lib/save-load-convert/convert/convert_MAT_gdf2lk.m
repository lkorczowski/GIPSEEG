function EEG= convert_MAT_gdf2lk(EEG_old)

EEG=EEG_old;
EEG.Trigger =EEG_old.Flash;
EEG.EpochClass=EEG_old.Y;
EEG.Channels=EEG_old.s;
EEG=removefields(EEG,{'Flash','Y','s'})
    try
    EEG.ElectrodesName=EEG_old.channels;
    end
end