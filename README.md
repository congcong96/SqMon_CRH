# SqMon_CRH
To get tmf and smf of stimulus files, run **batch_stimulus_to_tmf_smf**.  
**sm_calculate_CRH** uses binned spike trains to get CRH: 
the spike trains and the stimulus are binned at the same time resolution and the tmf and smf values are counted when a spike happens.  
There are two ways of calculating STRF:  
  **calculate_strf** uses .spr file of the stimulus and read blocks of stimulus for spike-triggered-average, thus the time resolution can be as fine as 0.05ms  
  **quick_calc_sta** uses matlab matrix product and requires loading the entire stimulus. So it is usually used for downsampled stimlus becasuse of the memory limit.




**ProbeRecord.mat** contains infromation about probe types and recording areas. Important fields: 
   * *monkey* - from which monkey the data is collected;   
   * *exp* - time stamp for the recording in the form yymmdd_hhmmss;   
   * *stim* - sound stimulus played during the recording;   
   * *ProbeName* - name of the probe used;   
   * *area* - recoridng made from core or belt, empty means no record available;   
   * *field* - fileds in core or belt area. We are less confident about the labels of fields within each area. So, be careful with your interpretation when using this.
