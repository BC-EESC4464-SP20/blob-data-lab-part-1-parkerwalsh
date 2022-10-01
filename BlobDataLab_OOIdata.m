%% 1. Explore and extract data from one year of OOI mooring data

clear
clc
load("matlab.mat")

filename = 'deployment0001_GP03FLMB.nc';
%1a. Use the function "ncdisp" to display information about the data contained in this file
%-->

ncdisp(filename);

%1b. Use the function "ncreadatt" to extract the latitude and longitude
%attributes of this dataset
%--> 
%-->

station_papa_lat_value = ncreadatt(filename,'/','lat');

station_papa_lon_value  = ncreadatt(filename,'/','lon');

%1c. Use the function "ncread" to extract the variables "time" and
%"ctdmo_seawater_temperature"
%-->
%-->

station_papa_time_values = ncread(filename,'time');
station_papa_seawater_temperature_values = ncread(filename,'ctdmo_seawater_temperature');
% pressure = ncread(filename,'ctdmo_seawater_pressure_qc_results');
% Extension option: Also extract the variable "pressure" (which, due to the
% increasing pressure underwater, tells us about depth - 1 dbar ~ 1 m
% depth). How deep in the water column was this sensor deployed?

%% 2. Converting the timestamp from the raw data to a format you can use
% Use the datenum function to convert the "time" variable you extracted
% into a MATLAB numerical timestamp (Hint: you will need to check the units
% of time from the netCDF file.)

% -->

station_papa_datenum_time_values = datenum(1900,1,1,0,0,station_papa_time_values);

% timemat = datenum(time);

% Checking your work: Use the "datestr" function to check that your
% converted times match the time range listed in the netCDF file's
% attributes for time coverage

station_papa_string_time_values = datestr(station_papa_datenum_time_values);

% 2b. Calculate the time resolution of the data (i.e. long from one
% measurement to the next) in minutes. Hint: the "diff" function will be
% helpful here.

res = mean(diff(station_papa_datenum_time_values));

%% 3. Make an initial exploration plot to investigate your data
% Make a plot of temperature vs. time, being sure to show each individual
% data point. What do you notice about the seasonal cycle? What about how
% the variability in the data changes over the year?
% Hint: Use the function "datetick" to make the time show up as
% human-readable dates rather than the MATLAB timestamp numbers

plot(station_papa_datenum_time_values, station_papa_seawater_temperature_values)
xlabel("Time")
ylabel("degC")
title("Seawater Temperature depth=30m at SMB Station Papa 2013-14")
datetick('x',"yyyy/mmm")

%% 4. Dealing with data variability: smoothing and choosing a variability cutoff
% 4a. Use the movmean function to calculate a 1-day (24 hour) moving mean
% to smooth the data. Hint: you will need to use the time period between
% measurements that you calculated in 2b to determine the correct window
% size to use in the calculation.

% -->

station_papa_seawater_temperature_movmean = movmean(station_papa_seawater_temperature_values,(1/res));

% 4b. Use the movstd function to calculate the 1-day moving standard
% deviation of the data.

station_papa_seawater_temperature_movstd = movstd(station_papa_seawater_temperature_values,(1/res));

%% 5. Honing your initial investigation plot
% Building on the initial plot you made in #3 above, now add:
%5a. A plot of the 1-day moving mean on the same plot as the original raw data
figure
subplot(2,1,1)
plot(station_papa_datenum_time_values, station_papa_seawater_temperature_values)
hold on
plot(station_papa_datenum_time_values, station_papa_seawater_temperature_movmean)
datetick('x','yyyy/mmm')
hold off
xlabel("Time")
ylabel("degC")
title("Raw data + 1-day movemean 2013-14")

%5b. A plot of the 1-day moving standard deviation, on a separate plot
%underneath, but with the same x-axis (hint: you can put two plots in the
%same figure by using "subplot" and you can specify

subplot(2,1,2)
plot(station_papa_datenum_time_values, station_papa_seawater_temperature_movstd)
datetick('x','yyyy/mmm')
xlabel("Time")
ylabel("Standard Deviation")
title("1-day moving standard deviation")

%% 6. Identifying data to exclude from analysis
% Based on the plot above, you can see that there are time periods when the
% data are highly variable - these are time periods when the raw data won't
% be suitable for use to use in studying the Blob.

%6a. Based on your inspection of the data, select a cutoff value for the
%1-day moving standard deviation beyond which you will exclude the data
%from your analysis. Note that you will need to justify this choice in the
%methods section of your writeup for this lab.
cutoff=2;

%6b. Find the indices of the data points that you are not excluding based
%on the cutoff chosen in 6a.
index = find(station_papa_seawater_temperature_movstd <= cutoff);

%6c. Update your figure from #5 to add the non-excluded data as a separate
%plotted set of points (i.e. in a new color) along with the other data you
%had already plotted.

figure
subplot(2,1,1)
plot(station_papa_datenum_time_values,station_papa_seawater_temperature_values)
hold on
plot(station_papa_datenum_time_values,station_papa_seawater_temperature_movmean,'k')
plot(station_papa_datenum_time_values(index),station_papa_seawater_temperature_movmean(index),'.r')
datetick('x','yyyy/mmm')
hold off
xlabel("Time")
ylabel('degC')
title("SST raw data + 2 std-cutoff smoothed data 2013-14")

subplot(2,1,2)
plot(station_papa_datenum_time_values,station_papa_seawater_temperature_movstd)
hold on
plot(station_papa_datenum_time_values(index),station_papa_seawater_temperature_movstd(index),'.r')
xlim([station_papa_datenum_time_values(1) station_papa_datenum_time_values(end)])
ylim([0 cutoff+0.5])
datetick('x','yyyy/mmm')
hold off
xlabel("Time")
ylabel("Standard Deviation")
title("Standard Deviation 2-std cutoff 2013-14")

%% 7. Apply the approach from steps 1-6 above to extract data from all OOI deployments in years 1-6
% You could do this by writing a for-loop or a function to adapt the code
% you wrote above to something you can apply across all 5 netCDF files
% (note that deployment 002 is missing)

concatenated_station_papa_seawater_temperatures = [];
concatenated_station_papa_datenum_times = [];


for i = [1,3:6]
loop_filename = ['deployment000' num2str(i) '_GP03FLMB.nc'];
loop_station_papa_time_values = ncread(loop_filename,'time');
loop_station_papa_seawater_temperature_values = ncread(loop_filename,'ctdmo_seawater_temperature');

loop_station_papa_datenum_time_values = datenum(1900,1,1,0,0,loop_station_papa_time_values);

res = mean(diff(loop_station_papa_datenum_time_values));


loop_station_papa_seawater_temperature_movmean = movmean(loop_station_papa_seawater_temperature_values,(1/res));

loop_station_papa_seawater_temperature_movstd = movstd(loop_station_papa_seawater_temperature_values,(1/res));

cutoff=2;

index = find(loop_station_papa_seawater_temperature_movstd <= cutoff);


figure
subplot(2,1,1)
plot(loop_station_papa_datenum_time_values,loop_station_papa_seawater_temperature_values)
hold on
plot(loop_station_papa_datenum_time_values,loop_station_papa_seawater_temperature_movmean,'k')
plot(loop_station_papa_datenum_time_values(index),loop_station_papa_seawater_temperature_movmean(index),'.r')
datetick
d_1 = datestr(loop_station_papa_datenum_time_values(1),'yyyy');
d_2 = datestr(loop_station_papa_datenum_time_values(end),'yyyy');
hold off
xlabel("Time")
ylabel('degC')
title("Raw SST + 1-day movemean at SFM B — OOI Station Papa array — " + d_1 + " - " + d_2)

subplot(2,1,2)
plot(loop_station_papa_datenum_time_values,loop_station_papa_seawater_temperature_movstd)
hold on
plot(loop_station_papa_datenum_time_values(index),loop_station_papa_seawater_temperature_movstd(index),'.r')
xlim([loop_station_papa_datenum_time_values(1) loop_station_papa_datenum_time_values(end)])
ylim([0 cutoff+0.5])
datetick
hold off
xlabel("Time")
ylabel("Standard Deviation")
title("Standard Deviation with "+string(cutoff)+"-std cutoff — " + d_1 + " - " + d_2)

if i == 1
concatenated_station_papa_datenum_times = [concatenated_station_papa_datenum_times; loop_station_papa_datenum_time_values(1:10:end);nan];
concatenated_station_papa_seawater_temperatures = [concatenated_station_papa_seawater_temperatures; loop_station_papa_seawater_temperature_movmean(1:10:end);nan]; 

else
concatenated_station_papa_seawater_temperatures = [concatenated_station_papa_seawater_temperatures; loop_station_papa_seawater_temperature_movmean(1:10:end)]; 
concatenated_station_papa_datenum_times = [concatenated_station_papa_datenum_times; loop_station_papa_datenum_time_values(1:10:end)];
end
end