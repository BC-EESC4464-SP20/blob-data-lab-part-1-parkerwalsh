%% 1. Explore and extract data from one year of OOI mooring data

clear
clc

filename = 'deployment0001_GP03FLMB.nc';
%1a. Use the function "ncdisp" to display information about the data contained in this file
%-->

ncdisp(filename);

%1b. Use the function "ncreadatt" to extract the latitude and longitude
%attributes of this dataset
%--> 
%-->

lat = ncreadatt(filename,'/','lat');

lon = ncreadatt(filename,'/','lon');

%1c. Use the function "ncread" to extract the variables "time" and
%"ctdmo_seawater_temperature"
%-->
%-->

time = ncread(filename,'time');
stemp = ncread(filename,'ctdmo_seawater_temperature');
% Extension option: Also extract the variable "pressure" (which, due to the
% increasing pressure underwater, tells us about depth - 1 dbar ~ 1 m
% depth). How deep in the water column was this sensor deployed?

%% 2. Converting the timestamp from the raw data to a format you can use
% Use the datenum function to convert the "time" variable you extracted
% into a MATLAB numerical timestamp (Hint: you will need to check the units
% of time from the netCDF file.)

% -->

timemat = datenum(1900,1,1,0,0,time);

% timemat = datenum(time);

% Checking your work: Use the "datestr" function to check that your
% converted times match the time range listed in the netCDF file's
% attributes for time coverage

timeconv = datestr(timemat);

% 2b. Calculate the time resolution of the data (i.e. long from one
% measurement to the next) in minutes. Hint: the "diff" function will be
% helpful here.

res = mean(diff(timemat));

%% 3. Make an initial exploration plot to investigate your data
% Make a plot of temperature vs. time, being sure to show each individual
% data point. What do you notice about the seasonal cycle? What about how
% the variability in the data changes over the year?
% Hint: Use the function "datetick" to make the time show up as
% human-readable dates rather than the MATLAB timestamp numbers

plot(timemat,stemp)
datetick

%% 4. Dealing with data variability: smoothing and choosing a variability cutoff
% 4a. Use the movmean function to calculate a 1-day (24 hour) moving mean
% to smooth the data. Hint: you will need to use the time period between
% measurements that you calculated in 2b to determine the correct window
% size to use in the calculation.

% -->

stempmm = movmean(stemp,(1/res));

% 4b. Use the movstd function to calculate the 1-day moving standard
% deviation of the data.

stempmstd = movstd(stemp,(1/res));

%% 5. Honing your initial investigation plot
% Building on the initial plot you made in #3 above, now add:
%5a. A plot of the 1-day moving mean on the same plot as the original raw data
figure
subplot(2,1,1)
plot(timemat,stemp)
hold on
plot(timemat,stempmm)
datetick
hold off

subplot(2,1,2)
plot(timemat,stempmstd)
datetick
%5b. A plot of the 1-day moving standard deviation, on a separate plot
%underneath, but with the same x-axis (hint: you can put two plots in the
%same figure by using "subplot" and you can specify

%% 6. Identifying data to exclude from analysis
% Based on the plot above, you can see that there are time periods when the
% data are highly variable - these are time periods when the raw data won't
% be suitable for use to use in studying the Blob.

%6a. Based on your inspection of the data, select a cutoff value for the
%1-day moving standard deviation beyond which you will exclude the data
%from your analysis. Note that you will need to justify this choice in the
%methods section of your writeup for this lab.

index = find(stempmstd <= 1);

%6b. Find the indices of the data points that you are not excluding based
%on the cutoff chosen in 6a.

%6c. Update your figure from #5 to add the non-excluded data as a separate
%plotted set of points (i.e. in a new color) along with the other data you
%had already plotted.

figure
subplot(2,1,1)
plot(timemat,stemp)
hold on
plot(timemat,stempmm,'k')
plot(timemat(index),stempmm(index),'.r')
datetick
hold off

subplot(2,1,2)
plot(timemat,stempmstd)
hold on
plot(timemat(index),stempmstd(index),'.r')
xlim([timemat(1) timemat(end)])
ylim([0 2.5])
datetick
hold off

%% 7. Apply the approach from steps 1-6 above to extract data from all OOI deployments in years 1-6
% You could do this by writing a for-loop or a function to adapt the code
% you wrote above to something you can apply across all 5 netCDF files
% (note that deployment 002 is missing)

for i = [1,3:6]
filename = ['deployment000' num2str(i) '_GP03FLMB.nc'];
time = ncread(filename,'time');
stemp = ncread(filename,'ctdmo_seawater_temperature');
timemat = datenum(1900,1,1,0,0,time);
res = mean(diff(timemat));
stempmm = movmean(stemp,(1/res));
stempmstd = movstd(stemp,(1/res));
index = find(stempmstd <= 1.5);

figure
subplot(2,1,1)
plot(timemat,stemp)
hold on
plot(timemat,stempmm,'k')
plot(timemat(index),stempmm(index),'.r')
datetick
hold off

subplot(2,1,2)
plot(timemat,stempmstd)
hold on
plot(timemat(index),stempmstd(index),'.r')
xlim([timemat(1) timemat(end)])
ylim([0 2.5])
datetick
hold off

end