load("matlab.mat")

%% Mapping Prep
imagescn(satellite_longitude, satellite_latitude, satellite_SST_anomaly(:,:,7)')
cmocean('thermal')

% index = find(mask ~= 1);

% sstAnom(index) = nan;
%% Mapping
figure

lat_limits = [min(satellite_latitude) max(satellite_latitude)];
lon_limits = [min(satellite_longitude) max(satellite_longitude)];

month_since_jan_2013 = 13;

axesm("MapProjection","mercator","MapLatLimit",lat_limits,"MapLonLimit",lon_limits);

pcolorm(satellite_latitude, satellite_longitude, satellite_SST_anomaly(:,:,month_since_jan_2013)');

geoshow('landareas.shp',"FaceColor","black");

caxis([-4 4]);
a = colorbar;
cmocean('balance','pivot');

a.Label.String="degC";
a.Label.Rotation=0;
a.Label.FontSize=12;
a.Label.Position = [2.75,0.5,0];

scatterm(50.5,-144.5,100,'k','o','LineWidth',2);
title("SST Anomalies in NE Pacific " + datestr(satellite_datenum_time_values(month_since_jan_2013),'mmm-yyyy'));
tightmap();
% grid = gridm("on");
% grid.Color = 'black';
% grid.LineWidth = 2.0;

%% Map and Line Subplots

figure;
subplot(1,2,1);

lat_limits = [min(satellite_latitude) max(satellite_latitude)];
lon_limits = [min(satellite_longitude) max(satellite_longitude)];

month_since_jan_2013 = 12;

axesm("MapProjection","mercator","MapLatLimit",lat_limits,"MapLonLimit",lon_limits);

pcolorm(satellite_latitude, satellite_longitude, satellite_SST_anomaly(:,:,month_since_jan_2013)');

geoshow('landareas.shp',"FaceColor","black");

caxis([-4 4]);
a = colorbar;
cmocean('balance','pivot');

a.Label.String="degC";
a.Label.Rotation=0;
a.Label.FontSize=12;
a.Label.Position = [2.75,0.5,0];

scatterm(50.5,-144.5,100,'k','o','LineWidth',2);
title("Satellite Seawater Temperature Anomalies in NE Pacific " + datestr(satellite_datenum_time_values(month_since_jan_2013),'mmm-yyyy'));
tightmap()

subplot(1,2,2);
plot(satellite_datenum_time_values, satellite_temperatures_at_station_papa);
hold on
plot(concatenated_station_papa_datenum_times, station_papa_seawater_temperature_anomaly,Color=[0.8500 0.3250 0.0980 0.3])

datetick('x','yyyy/mmm')
ylabel("degC")
xlabel("time")
title("Comparison of buoy and satellite temperature anomaly at depth = 30m")
legend(["Satellite data" "Buoy data"],"Box","off","AutoUpdate","off")

xlim([satellite_datenum_time_values(1) satellite_datenum_time_values(end)])

grid("on")
grid("minor")

yline(0)

xline(satellite_datenum_time_values(month_since_jan_2013))

% xticks(satellite_datenum_time_values(1:10:end))


% zero_values = zeros(size(satellite_datenum_time_values));

% plot(satellite_datenum_time_values, zero_values,'k','LineWidth',0.7)

camroll(90)

%% GIF making

figure;

axesm("MapProjection","pcarree","MapLatLimit",lat_limits,"MapLonLimit",lon_limits)

pcolorm(lat,lon,sstAnom(:,:,1)')
% % coast = load("coastlines.mat");
% geoshow(coast.coastlat,coast.coastlon)
geoshow('landareas.shp',"FaceColor","black")
% geoshow(lat,lon,sstAnom(:,:,52),cmocean('balance','pivot')
caxis([-4 4])
a = colorbar;
a.Label.String="degC"
a.Label.Rotation=0
a.Label.FontSize=12
a.Label.Position = [2.75,0.5,0]
cmocean('balance','pivot')

scatterm(50.5,-144.5,100,'yellow','filled','o')
title("SST Anomalies in NE Pacific " + datestr(time(1),'mmm-yyyy'))


gif('myGif.gif','DelayTime',1/6)

for i = 2:length(time)
    
    axesm("MapProjection","pcarree","MapLatLimit",latlim,"MapLonLimit",lonlim)

    pcolorm(lat,lon,sstAnom(:,:,i)')
    
% % coast = load("coastlines.mat");
% geoshow(coast.coastlat,coast.coastlon)
    geoshow('landareas.shp',"FaceColor","black")
% geoshow(lat,lon,sstAnom(:,:,52),cmocean('balance','pivot')
    caxis([-4 4])
    a = colorbar;
    a.Label.String="degC"
    a.Label.Rotation=0
    a.Label.FontSize=12
    a.Label.Position = [2.75,0.5,0]
    cmocean('balance','pivot')

    scatterm(50.5,-144.5,100,'yellow','filled','o')
    title("SST Anomalies in NE Pacific " + datestr(time(i),'mmm-yyyy'))
    gif
end
% figure; 
% axesm miller
% h = geoshow('landareas.shp');

%% Gif 2.0
figure;
subplot(1,2,1);

lat_limits = [min(satellite_latitude) max(satellite_latitude)];
lon_limits = [min(satellite_longitude) max(satellite_longitude)];

month_since_jan_2013 = 1;

axesm("MapProjection","mollweid","MapLatLimit",lat_limits,"MapLonLimit",lon_limits);

pcolorm(satellite_latitude, satellite_longitude, satellite_SST_anomaly(:,:,month_since_jan_2013)');

geoshow('landareas.shp',"FaceColor","black");

caxis([-4 4]);
a = colorbar;
cmocean('balance','pivot');

a.Label.String="degC";
a.Label.Rotation=0;
a.Label.FontSize=12;
a.Label.Position = [2.75,0.5,0];

scatterm(50.5,-144.5,100,'k','o','LineWidth',2);
title("Satellite Seawater Temperature Anomalies in NE Pacific " + datestr(satellite_datenum_time_values(month_since_jan_2013),'mmm-yyyy'));
tightmap()

subplot(1,2,2);
plot(satellite_datenum_time_values, satellite_temperatures_at_station_papa);
hold on
plot(concatenated_station_papa_datenum_times, station_papa_seawater_temperature_anomaly,Color=[0.8500 0.3250 0.0980 0.3])

datetick('x','yyyy/mmm')
ylabel("degC")
xlabel("time")
title("Comparison of buoy and satellite temperature anomaly at depth = 30m")
legend(["Satellite data" "Buoy data"],"Box","off","AutoUpdate","off")

xlim([satellite_datenum_time_values(1) satellite_datenum_time_values(end)])

grid("on")
grid("minor")

yline(0)

xline(satellite_datenum_time_values(month_since_jan_2013))

% xticks(satellite_datenum_time_values(1:10:end))


% zero_values = zeros(size(satellite_datenum_time_values));

% plot(satellite_datenum_time_values, zero_values,'k','LineWidth',0.7)

camroll(90)

gif('myGif1_4sat_moll.gif','DelayTime',1/4)
for month_since_jan_2013 = 2:length(satellite_datenum_time_values)
    clf
    subplot(1,2,1);

    lat_limits = [min(satellite_latitude) max(satellite_latitude)];
    lon_limits = [min(satellite_longitude) max(satellite_longitude)];
    
    % month_since_jan_2013 = 1;
    
    axesm("MapProjection","mollweid","MapLatLimit",lat_limits,"MapLonLimit",lon_limits);
    
    pcolorm(satellite_latitude, satellite_longitude, satellite_SST_anomaly(:,:,month_since_jan_2013)');
    
    geoshow('landareas.shp',"FaceColor","black");
    
    caxis([-4 4]);
    a = colorbar;
    cmocean('balance','pivot');
    
    a.Label.String="degC";
    a.Label.Rotation=0;
    a.Label.FontSize=12;
    a.Label.Position = [2.75,0.5,0];
    
    scatterm(50.5,-144.5,100,'k','o','LineWidth',2);
    title("Satellite Seawater Temperature Anomalies in NE Pacific " + datestr(satellite_datenum_time_values(month_since_jan_2013),'mmm-yyyy'));
    tightmap()
    
    subplot(1,2,2);
    plot(satellite_datenum_time_values, satellite_temperatures_at_station_papa);
    hold on
    plot(concatenated_station_papa_datenum_times, station_papa_seawater_temperature_anomaly,Color=[0.8500 0.3250 0.0980 0.3])
    
    datetick('x','yyyy/mmm')
    ylabel("degC")
    xlabel("time")
    title("Comparison of buoy and satellite temperature anomaly at depth = 30m")
    legend(["Satellite data" "Buoy data"],"Box","off","AutoUpdate","off")
    
    xlim([satellite_datenum_time_values(1) satellite_datenum_time_values(end)])
    
    grid("on")
    grid("minor")
    
    yline(0)
    
    xline(satellite_datenum_time_values(month_since_jan_2013))
    
    % xticks(satellite_datenum_time_values(1:10:end))
    
    
    % zero_values = zeros(size(satellite_datenum_time_values));
    
    % plot(satellite_datenum_time_values, zero_values,'k','LineWidth',0.7)
    
    camroll(90)
    gif
end