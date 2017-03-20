#read in the raw data for PHICS and return the blip data for the infants (basically all positive swabs before infection date)

#Files that this saves out:
## 1. CMVblipdata, CMVblipdata_pos -> cmv_pos_data_altdef.RData: The old infection definition
## 2. all_blips -> all_pos_data.RData: The pre-infection period of blips with correct right censoring
## 3. CMVblipdata, blip_only_data_CMV -> cmv_blip_data.RData: Pre-infection blip period and then blip-only (no negative swabs) -CMV only
## 4. all_blips, blip_only_data -> all_blip_data.RData: Pre-infection blip period and then blip-only (no negative swabs) - all viruses
## 5. inf_data, swab_size_data -> infected_blip_data.RData: CMV data before and after infection time for comparison
## 6. blip_data_duration, blip_data_duration_complete -> duration_data.RData: Blips with duration information, complete: excluded censored

library(dplyr)

virusMeltedDataDemoAllInfant = readr::read_csv("data/raw_data.csv")

infantCMVdata = subset(virusMeltedDataDemoAllInfant, Virus == "ORL_CMV") %>%
    group_by(PatientID2) %>% mutate(last_date = unique(infantInfDate))

infantCMVdata$keep_blips = with(infantCMVdata, ifelse(infantInfection == 1, (times < last_date), TRUE))

CMVblipdata = select(subset(infantCMVdata, keep_blips), -keep_blips) %>% ungroup()
CMVblipdata$pos = CMVblipdata$count > 0
CMVblipdata = CMVblipdata %>% group_by(PatientID2) %>%
  mutate(time_diff = c(NA, diff(days)),
         consecutive_swab = lag(pos) & pos)

#----------------- make all transient infection data -----------------------------
other_virus_blips = subset(virusMeltedDataDemoAllInfant, Virus != "ORL_CMV")
other_virus_blips$keep = with(other_virus_blips, ifelse(infantInfection == 1, (times < infantInfDate), TRUE))
all_blips_temp = select(subset(other_virus_blips, keep), -keep) %>% group_by(PatientID2, Virus) %>%
  mutate(pos = count > 0,
         time_diff = c(NA, diff(days)),
         consecutive_swab = lag(pos) & pos) %>%
  bind_rows(CMVblipdata) %>% ungroup()
all_blips_temp$Virus = matrix(unlist(strsplit(as.character(all_blips_temp$Virus), '_')), ncol=2, byrow=T) [, 2]


###### Need to positive swabs at the end of the study that are not followed by a negative swab

all_blips = plyr::ldply(unique(all_blips_temp$Virus), function(virus){
  sub_virus_data = subset(all_blips_temp, Virus == virus)
  plyr::ldply(unique(sub_virus_data$PatientID2), function(pid){
    sub_data = subset(sub_virus_data, PatientID2 == pid)
    sub_data = sub_data %>% mutate(check = (times == max(times) & pos))
    while(any(sub_data$check)) sub_data = filter(sub_data, !check) %>% mutate(check = (times == max(times) & pos))
    select(sub_data, -check)
  })
})

save(all_blips, file = "data/all_pos_data.RData")

#------------------------- create data with only the transient infections and they are counted ------------------#
#crawls through the data and isolated consecutive blips
define_blip = function(blip_no, blip_data, prev_day, next_day, debug = F){
  if(debug) browser()
  data.frame(
    blip_no = blip_no,
    count = blip_data$count,
    days = blip_data$days,
    prev_day = c(prev_day, head(blip_data$days, -1)),
    next_day = c(tail(blip_data$days, -1), next_day)
  )
}

create_blips = function(data_in, debug = F){
  if(debug) browser()
  viral_ts = select(data_in, days, count) %>% arrange(days)
  total_swabs = dim(viral_ts)[1]

  if(total_swabs == 0) {print(paste(data_in$PatientID2[1], data_in$Virus[1], "has length 0")); return(NULL)} #no blips (weird case)
  if(total_swabs == 1){
    if(viral_ts$count[1] > 0) return(define_blip(1, viral_ts, NA, NA)) else return(NULL)
  }

  out_blip_data = NULL
  blip_no = 1
  prev_day = NA

  is_blip = viral_ts$count[1] > 0
  if(is_blip) start_index = 1

  for(i in 2:total_swabs){
    if(is_blip){
      if(viral_ts$count[i] == 0){ #blip ends
        out_blip_data = bind_rows(out_blip_data,
                                  define_blip(blip_no, viral_ts[start_index:(i - 1), ],
                                              prev_day, viral_ts$days[i])
        )
        blip_no = blip_no + 1
        is_blip = F
      }
    } else{
      if(viral_ts$count[i] > 0){
        is_blip = T
        prev_day = viral_ts$days[i - 1]
        start_index = i
      }
    }
  }
  #if final censored blip
  if(is_blip) out_blip_data = bind_rows(out_blip_data,
                                        define_blip(blip_no, viral_ts[start_index:total_swabs, ],
                                                    prev_day, NA)
  )

  out_blip_data
}

blip_only_data =  plyr::ldply(unique(all_blips$Virus), function(v){
  if(v %in% c("HHV8", "HSV")) return(NULL)
  print(v)
  temp_data = subset(all_blips, Virus == v)
  plyr::ldply(unique(temp_data$PatientID2), function(pid){
    out = create_blips(subset(temp_data, PatientID2 == pid))
    if(!is.null(out)) {out$PatientID2 = pid; out$Virus = v; out$infection = subset(temp_data, PatientID2 == pid)$infantInfection[1]}
    out
    })
})

#should have 297 - (16 + 70 from removed HHV8 and HSV) = 211 blips (it checks)
blip_only_data %>% group_by(PatientID2, Virus, blip_no) %>% summarize() %>% count(Virus) %>% mutate(sum(n))

#load the Elizabeth check portion to double check
#my_included %>% group_by(Virus) %>% summarize(total = n())
#EK_included %>% group_by(Virus) %>% summarize(total = n())

blip_only_data_CMV = subset(blip_only_data, Virus == "CMV")
blip_only_data_CMV %>% group_by(PatientID2, blip_no) %>% summarize() %>% count()

save(CMVblipdata, blip_only_data_CMV, file = "data/cmv_blip_data.RData")
save(all_blips, blip_only_data, file = "data/all_blip_data.RData")


#----------- Make data to compare infected shedding data to blip shedding data in CMV ------------

CMVPrimaryEpisodes = readr::read_csv("data/CMV_primary_infection_data.csv")


excluded_IDs = setdiff(unique(subset(infantCMVdata, infantInfection == 1)$PatientID2), unique(CMVPrimaryEpisodes$PatientID2))
excluded_ID_data = subset(infantCMVdata, PatientID2 %in% excluded_IDs & times >= infantInfDate)

inf_data = bind_rows(CMVPrimaryEpisodes, subset(excluded_ID_data, PatientID2 != "02537-P"))

swab_size_data = bind_rows(blip_only_data_CMV %>% group_by(PatientID2) %>% summarize(type = "Blip", max_size = max(count, na.rm = T)),
                           inf_data %>% group_by(PatientID2) %>% summarize(type = "Infection", max_size = max(count, na.rm = T)))


save(inf_data, swab_size_data, file = "data/infected_blip_data.RData")

#-------  Duration data -----#
blip_data_duration = blip_only_data_CMV %>% group_by(PatientID2, blip_no) %>%
  summarize(
    min_day_diff = min(c(diff(days), min(days) - min(prev_day), max(next_day) - max(days)), na.rm = T),
    max_day_diff = max(c(diff(days), min(days) - min(prev_day), max(next_day) - max(days)), na.rm = T),
    start_day_est = min(prev_day) + (min(days) - min(prev_day))/2,
    end_day_est =  max(days) + (max(next_day) - max(days))/2,
    duration = end_day_est - start_day_est,
    total_blips = n())
blip_data_duration_complete = subset(blip_data_duration, min_day_diff >= 5 & max_day_diff <= 9)

blip_data_duration_excluded = subset(blip_data_duration, min_day_diff < 5 | max_day_diff > 9)
hist(blip_data_duration_excluded$total_blips)

sum(blip_data_duration_excluded$total_blips)

blip_data_duration_excluded %>% ungroup() %>% dplyr::mutate(total = n()) %>% group_by(total_blips) %>% summarize(n()/total[1])

save(blip_data_duration, blip_data_duration_complete, file = "data/duration_data.RData")


