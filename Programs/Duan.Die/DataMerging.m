date_beg = datetime(2020,7,26);
date_end = datetime(2020,10,27);
temp_date = date_beg-1;
fname_mag = 'solo_l2_mag-rtn-normal_';
fname_pas = 'solo_l2_swa-pas-grnd-mom_';
data_dir = 'D:\SolOData\';
epoch_tot = [];
B_tot = [];
np_tot = [];
Tp_tot = [];
r_tot = [];
Vp_tot = [];
posfile = 'D:\SolOData\solo_helio1day_position_20200211_v01.cdf';
pos_epoch = double(spdfcdfread(posfile,'Variables',{'Epoch'}));
r_all = double(spdfcdfread(posfile,'Variables',{'RAD_AU'}));
while temp_date<date_end
    temp_date = temp_date+1;
    temp_date_str = datestr(temp_date,'yyyymmdd');
    tempf_mag = dir([data_dir fname_mag temp_date_str '_v0*.cdf']);
    tempf_pas = dir([data_dir fname_pas temp_date_str '_v0*.cdf']);
    if isempty(tempf_pas)
        continue
    end
    if isempty(tempf_mag)
        continue
    end
    momsdir = [data_dir tempf_pas.name];
    magdir = [data_dir tempf_mag.name];
    epochmag = spdfcdfread(magdir,'Variables',{'EPOCH'});
B_rtn = double(spdfcdfread(magdir,'Variables',{'B_RTN'}));
epochpas = spdfcdfread(momsdir,'Variables',{'Epoch'});
N_rtn = double(spdfcdfread(momsdir,'Variables',{'N'}));
T_rtn = double(spdfcdfread(momsdir,'Variables',{'T'}));
V_rtn = double(spdfcdfread(momsdir,'Variables',{'V_RTN'}));
B_temp = interp1(epochmag,B_rtn,epochpas);
epoch_tot = [epoch_tot;epochpas];
B_tot = [B_tot;B_temp];
np_tot = [np_tot;N_rtn];
Tp_tot = [Tp_tot;T_rtn];
Vp_tot = [Vp_tot;V_rtn];
end
r_tot = interp1(pos_epoch,r_all,epoch_tot);
save('Jingyu.mat','epoch_tot','B_tot','Vp_tot','np_tot','Tp_tot','r_tot')