function [ state, FRU_state, Ltype_state, RyR_state, Ito2_state ] = initialize(ic_states_file, ic_FRU_file, ic_LCh_file,ic_RyR_file,ic_Ito2_file)

% Read in initial conditions from specified text files 
 global Nclefts_FRU NRyRs_per_cleft Nstates_FRU Nindepstates_LType
 global NFRU_sim N N_ss
% Define the indices
index_V         = 1;
index_mNa       = 2;
index_hNa       = 3;
index_jNa       = 4;
index_Nai       = 5;
index_Ki        = 6;
index_Cai       = 7;
index_CaNSR     = 8;
index_xKs       = 9;
index_LTRPNCa   = 10;
index_HTRPNCa   = 11;
index_C0Kv43    = 12;
index_C1Kv43    = 13;
index_C2Kv43    = 14;
index_C3Kv43    = 15;
index_OKv43     = 16;
index_CI0Kv43   = 17;
index_CI1Kv43   = 18;
index_CI2Kv43   = 19;
index_CI3Kv43   = 20;
index_OIKv43    = 21;
index_C0Kv14    = 22;
index_C1Kv14    = 23;
index_C2Kv14    = 24;
index_C3Kv14    = 25;
index_OKv14     = 26;
index_CI0Kv14   = 27;
index_CI1Kv14   = 28;
index_CI2Kv14   = 29;
index_CI3Kv14   = 30;
index_OIKv14    = 31;
index_CaTOT     = 32;
index_C1Herg    = 33;
index_C2Herg    = 34;
index_C3Herg    = 35;
index_OHerg     = 36;
index_IHerg     = 37;

state       = zeros(N,1);
FRU_state   = zeros(N_ss, Nstates_FRU);
Ltype_state = zeros(N_ss, Nclefts_FRU,Nindepstates_LType);
RyR_state   = zeros(NFRU_sim, Nclefts_FRU,NRyRs_per_cleft);
Ito2_state  = zeros(N_ss,Nclefts_FRU);

% Initialize the states
fileID = fopen(ic_states_file,'r');
C_data  = textscan(fileID,'%s');
state(index_V)        = str2double(C_data{1}(index_V));
state(index_mNa)      = str2double(C_data{1}(index_mNa));
state(index_hNa)      = str2double(C_data{1}(index_hNa));
state(index_jNa)      = str2double(C_data{1}(index_jNa));
state(index_Nai)      = str2double(C_data{1}(index_Nai));
state(index_Ki)       = str2double(C_data{1}(index_Ki));
state(index_Cai)      = str2double(C_data{1}(index_Cai));
state(index_CaNSR)    = str2double(C_data{1}(index_CaNSR));
state(index_LTRPNCa)  = str2double(C_data{1}(index_LTRPNCa));
state(index_HTRPNCa)  = str2double(C_data{1}(index_HTRPNCa));
state(index_xKs)      = str2double(C_data{1}(index_xKs));
state(index_C1Herg)   = str2double(C_data{1}(index_C1Herg));
state(index_C2Herg)   = str2double(C_data{1}(index_C2Herg));
state(index_C3Herg)   = str2double(C_data{1}(index_C3Herg));
state(index_OHerg)    = str2double(C_data{1}(index_OHerg));
state(index_IHerg)    = str2double(C_data{1}(index_IHerg));
state(index_C0Kv43)   = str2double(C_data{1}(index_C0Kv43));
state(index_C1Kv43)   = str2double(C_data{1}(index_C1Kv43));
state(index_C2Kv43)   = str2double(C_data{1}(index_C2Kv43));
state(index_C3Kv43)   = str2double(C_data{1}(index_C3Kv43));
state(index_OKv43)    = str2double(C_data{1}(index_OKv43));
state(index_CI0Kv43)  = str2double(C_data{1}(index_CI0Kv43));
state(index_CI1Kv43)  = str2double(C_data{1}(index_CI1Kv43));
state(index_CI2Kv43)  = str2double(C_data{1}(index_CI2Kv43));
state(index_CI3Kv43)  = str2double(C_data{1}(index_CI3Kv43));
state(index_OIKv43)  = 1.0- sum(state(index_OIKv43));
state(index_C0Kv14)  = str2double(C_data{1}(index_C0Kv14));
state(index_C1Kv14)  = str2double(C_data{1}(index_C1Kv14));
state(index_C2Kv14)  = str2double(C_data{1}(index_C2Kv14));
state(index_C3Kv14)   = str2double(C_data{1}(index_C3Kv14));
state(index_OKv14)   = str2double(C_data{1}(index_OKv14));
state(index_CI0Kv14)  = str2double(C_data{1}(index_CI0Kv14));
state(index_CI1Kv14)  = str2double(C_data{1}(index_CI1Kv14));
state(index_CI2Kv14)  = str2double(C_data{1}(index_CI2Kv14));
state(index_CI3Kv14) = str2double(C_data{1}(index_CI3Kv14));
state(index_OIKv14)  =  1.0- sum(state(index_OIKv14));
state(index_CaTOT)  =   str2double(C_data{1}(index_CaTOT));

state = state';

% Initialize the FRU concentrations

fileID = fopen(ic_FRU_file, 'r');
C_data1  = textscan(fileID,'%s %s %s %s %s');
IC_length = length(str2double(C_data1{1}(2:end)));
rewind_times = ceil(N_ss/IC_length);
i = 1;
while (i <= rewind_times)
    start_ind = (i-1)*250 + 1;
    end_ind = min(i*250, N_ss);
    FRU_state(start_ind:end_ind,1) = str2double(C_data1{1}(2:(2+ end_ind - start_ind)));
    FRU_state(start_ind:end_ind,2)  = str2double(C_data1{2}(2:(2+ end_ind - start_ind)));
    FRU_state(start_ind:end_ind,3)  = str2double(C_data1{3}(2:(2+ end_ind - start_ind)));
    FRU_state(start_ind:end_ind,4)  = str2double(C_data1{4}(2:(2+ end_ind - start_ind)));
    FRU_state(start_ind:end_ind,5)  = str2double(C_data1{5}(2:(2+ end_ind - start_ind)));
    i = i+1;
end

fileID = fopen(ic_LCh_file);
C_data2  = textscan(fileID,'%s %s');
IC_length = IC_length/4;
rewind_times = ceil(N_ss/IC_length);
L = str2double(C_data2{1}(2:end));
Y = str2double(C_data2{2}(2:end));
r = 1;
while(r<= rewind_times)
  i=1;
  FRU_start = (r-1)*250 + 1;
  FRU_end = min((r*250),N_ss);
  for iFRU=FRU_start:FRU_end
        for icleft = 1:Nclefts_FRU
               Ltype_state(iFRU,icleft,1) = L(i);
               Ltype_state(iFRU,icleft,2) = Y(i);
               i = i+ 1;
        end
  end
  r =  r+ 1;
end

fileID = fopen(ic_RyR_file);
C_data3  = textscan(fileID,'%s');
M = str2double(C_data3{1}(2:end));
IC_length = length(str2double(C_data2{1}(2:end)));
IC_length = IC_length/(Nclefts_FRU*NRyRs_per_cleft);
rewind_times = ceil(NFRU_sim/IC_length);
r = 1;
while(r< rewind_times)
    i = 1;
    FRU_start = (r-1)*250 + 1;
    FRU_end = min((r*250),NFRU_sim);
    for iFRU=FRU_start:FRU_end
        for icleft=1:Nclefts_FRU
            for iRyR=1:NRyRs_per_cleft
                RyR_state(iFRU,icleft,iRyR) = M(i);
                i = i+ 1;
            end
         end
    end
    r = r+1;
end

fileID = fopen(ic_Ito2_file );
C_data4  = textscan(fileID,'%s');
IC_length = length(str2double(C_data2{1}(2:end)));
IC_length = IC_length/Nclefts_FRU;
rewind = ceil(N_ss/ IC_length);
r =1;
while(r<= rewind)
start = (r-1)*IC_length + 1;
end_ind = min(r*IC_length, N_ss);
Ito2_state(start:end_ind, 1) = str2double(C_data4{1}(2:(2+end_ind - start)));
Ito2_state(start:end_ind, 2) = str2double(C_data4{1}(252:(252+end_ind - start)));
Ito2_state(start:end_ind, 3) = str2double(C_data4{1}(502:(502+end_ind - start)));
Ito2_state(start:end_ind, 4) = str2double(C_data4{1}(752:(752+end_ind - start)));
r= r+1;
end

end

