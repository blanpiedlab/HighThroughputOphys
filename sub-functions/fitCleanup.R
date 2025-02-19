#fitCleanup.R




ROIpeaks_ThatPass<- unique(peaksThatPass$ROIpeaks)
fitsCleaned<- sync_async_fitsLabeled_stimKeyed %>% dplyr::filter(ROIpeaks %in% ROIpeaks_ThatPass)
