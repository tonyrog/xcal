%%% @author Tony Rogvall <tony@rogvall.se>
%%% @copyright (C) 2020, Tony Rogvall
%%% @doc
%%%    Calculate sunset & sunrise based on location 
%%%    converted from: https://github.com/kelvins/sunrisesunset
%%% @end
%%% Created : 11 Apr 2020 by Tony Rogvall <tony@rogvall.se>

-module(sunset_sunrise).

-export([get_time/0]).
-export([get_time/1]).
-export([get_time/4]).
-export([benchmark/0]).
-export([test/0, test/5]).

-compile(inline).
-compile({inline_size,100}).


-define(SECONDS_PER_DAY, (24*60*60)).

-type radians() :: float().
-type degree() :: float().
-type minute() :: float().
-type second() :: float().

%% Convert radians to degrees
-define(rad2deg(Radians), ((Radians)*(180.0 / math:pi()))).
-define(deg2rad(Degrees), ((Degrees)*(math:pi() / 180.0))).

sqr(X) -> X*X.

%% bench: 4.25  (2020-04-11 15:42)
%% bench: 4.98  (2020-04-11 16:14)
%% bench: 8.86  (2020-04-11 23:10)
%% bench: 10.35 (2020-04-12 00:03)
%% bench: 10.60 (2020-04-12 00:27)
%% bench: 11.16 (2020-04-12 01:26)

%% sunrisesunset is used to calculate the apparent sunrise and sunset 
%% based on the latitude, longitude, UTC offset and date.
%% All calculations (formulas) were extracted from the 
%% Solar Calculation Details of the Earth System Research Laboratory:
%% https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html
%%

%% get_time function is responsible for calculate the apparent 
%% Sunrise and Sunset times.
%% If some parameter is wrong it will return an error.

get_time() ->
    get_time(date()).
    
get_time(Date) ->
    {Latitude,Longitude,UtcOffset} = eplanet:location(),
    get_time(Latitude,Longitude,UtcOffset,Date).

-spec get_time(Latitude::number(), Longitude::number(), UtcOffset::number(), Date::calendar:date()) ->
		      {error, Message::string()} | 
		      {ok, {Sunrise::calendar:time(),Sunset::calendar:time()}}.

get_time(Latitude, Longitude, UtcOffset, Date={Year,_Month,_Day}) when
      (Latitude >= -90),(Latitude =< 90),
      (Longitude >= -180), (Longitude =< 180),
      (UtcOffset >= -12), (UtcOffset =< 14),
      (Year >= 1900), (Year =< 2200) ->
    case calendar:valid_date(Date) of
	false -> {error, invalid_date};
	true ->
	    get_time_(Latitude, Longitude, UtcOffset, Date)
    end.

%% I guess second accuracy is a bit over kill just to find
%% out when sun is settings, bit who knows?
get_time_(Latitude, Longitude, UtcOffset, Date) ->
    NumDays = days_diff({1899,12,30}, Date),
    {SunriseSeconds,SunsetSeconds} = 
	get_time__(Latitude, Longitude, UtcOffset, NumDays,
		   0, ?SECONDS_PER_DAY, 
		   0, ?SECONDS_PER_DAY, 0, ?SECONDS_PER_DAY),
    Sunrise = calendar:seconds_to_time(SunriseSeconds),
    Sunset = calendar:seconds_to_time(SunsetSeconds),
    {ok,{Sunrise, Sunset}}.

get_time__(Latitude, Longitude, UtcOffset, NumDays, 
	   Second, SecondMax,
	   Sunrise, SunriseMin,
	   Sunset, SunsetMin) when Second < SecondMax ->
    SecondNorm = Second / (?SECONDS_PER_DAY-1),
    JulianDay = calcJulianDay(NumDays, SecondNorm, UtcOffset),
    JulianCentury = calcJulianCentury(JulianDay),
    GeomMeanLongSun_r = calcGeomMeanLongSun_r(JulianCentury),
    GeomMeanAnomSun_r = calcGeomMeanAnomSun_r(JulianCentury),
    EccentEarthOrbit = calcEccentEarthOrbit(JulianCentury),
    SunEqCtr = calcSunEqCtr(JulianCentury, GeomMeanAnomSun_r),
    SunTrueLong = calcSunTrueLong(SunEqCtr, GeomMeanLongSun_r),
    SunAppLong_r = calcSunAppLong_r(SunTrueLong, JulianCentury),
    MeanObliqEcliptic = calcMeanObliqEcliptic(JulianCentury),
    ObliqCorr_r = calcObliqCorr_r(MeanObliqEcliptic, JulianCentury),
    SunDeclination_r = calcSunDeclination_r(ObliqCorr_r, SunAppLong_r),
    MultiFactor = sqr(math:tan(ObliqCorr_r/2.0)),
    EquationOfTime = calcEquationOfTime(MultiFactor, GeomMeanLongSun_r, 
					EccentEarthOrbit, GeomMeanAnomSun_r),
    HaSunrise = calcHaSunrise(Latitude, SunDeclination_r),
    SolarNoon = calcSolarNoon(Longitude, EquationOfTime, UtcOffset),

    Second0 = SecondNorm * ?SECONDS_PER_DAY,
    SunriseSecond = round(HaSunrise*4.0*60.0),
    Sunrise0   = SolarNoon - SunriseSecond - Second0,
    Sunset0    = SolarNoon + SunriseSecond - Second0,
    SunriseAbs = abs(Sunrise0),
    SunsetAbs  = abs(Sunset0),
    {Sunrise1,SunriseMin1} =
	if SunriseAbs < SunriseMin -> {Second, SunriseAbs};
	   true -> {Sunrise, SunriseMin}
	end,
    {Sunset1,SunsetMin1} =
	if SunsetAbs < SunsetMin -> {Second, SunsetAbs};
	   true -> {Sunset, SunsetMin}
	end,
    get_time__(Latitude, Longitude, UtcOffset, NumDays, 
	       Second+1, SecondMax,
	       Sunrise1, SunriseMin1,
	       Sunset1, SunsetMin1);
get_time__(_Latitude, _Longitude, _UtcOffset, _NumDays, 
	   Second, Second,
	   Sunrise, _SunriseMin,
	   Sunset, _SunsetMin) ->
    {Sunrise, Sunset}.

%%
%% Calculate Julian Day based on the formula:
%%   nDays+2415018.5+secondsNorm-UtcOffset/24
%% numDays - The number of days calculated in the calculate function
%% secondsNorm - Seconds normalized calculated by the createSecondsNormalized function
%% utcOffset - UTC offset defined by the user

calcJulianDay(NumDays, SecondsNorm, UtcOffset) ->
    NumDays + 2415018.5 + SecondsNorm - UtcOffset/24.0.

%% Calculate the Julian Century based on the formula: (julianDay - 2451545.0) / 36525.0
%% julianDay - Julian day vector calculated by the calcJulianDay function

calcJulianCentury(JulianDay) ->
    ((JulianDay - 2451545.0) / 36525.0).

%% Calculate the Geom Mean Long Sun in degrees based on the formula: 280.46646 + julianCentury * (36000.76983 + julianCentury * 0.0003032)
%% julianCentury - Julian century calculated by the calcJulianCentury function

-spec calcGeomMeanLongSun_r(JulianCentury::float()) -> radians().

calcGeomMeanLongSun_r(JulianCentury) ->
    ?deg2rad(math:fmod(280.46646 + JulianCentury*(36000.76983+JulianCentury*0.0003032),360.0)).

%% Calculate the Geom Mean Anom Sun in degrees based on the formula: 357.52911 + julianCentury * (35999.05029 - 0.0001537 * julianCentury)
%% julianCentury - Julian century calculated by the calcJulianCentury function

calcGeomMeanAnomSun_r(JulianCentury) ->
    ?deg2rad(357.52911 + JulianCentury*(35999.05029-0.0001537*JulianCentury)).

%% Calculate the Eccent Earth Orbit based on the formula: 0.016708634 - julianCentury * (0.000042037 + 0.0000001267 * julianCentury)
%% julianCentury - Julian century calculated by the calcJulianCentury function

calcEccentEarthOrbit(JulianCentury) ->
    (0.016708634 - JulianCentury*(0.000042037+0.0000001267*JulianCentury)).

%% Calculate the Sun Eq Ctr based on the formula: sin(deg2rad(geomMeanAnomSun))*(1.914602-julianCentury*(0.004817+0.000014*julianCentury))+sin(deg2rad(2.0*geomMeanAnomSun))*(0.019993-0.000101*julianCentury)+sin(deg2rad(3.0*geomMeanAnomSun))*0.000289;

-spec calcSunEqCtr(JulianCentury::float(),GeomMeanAnomSun_r::radians()) ->
			  float().

calcSunEqCtr(JulianCentury,GeomMeanAnomSun_r) ->
    math:sin(GeomMeanAnomSun_r)*(1.914602-JulianCentury*(0.004817+0.000014*JulianCentury)) +
    math:sin(2.0*GeomMeanAnomSun_r)*(0.019993-0.000101*JulianCentury) +
    math:sin(3.0*GeomMeanAnomSun_r)*0.000289.


%% Calculate the Sun True Long in degrees based on the formula: sunEqCtr + geomMeanLongSun
%% sunEqCtr - Sun Eq Ctr calculated by the calcSunEqCtr function
%% geomMeanLongSun - Geom Mean Long Sun calculated by the calcGeomMeanLongSun function

calcSunTrueLong(SunEqCtr, GeomMeanLongSun_r) ->
    (SunEqCtr + ?rad2deg(GeomMeanLongSun_r)).


%% Calculate the Sun App Long in degrees based on the formula: sunTrueLong-0.00569-0.00478*sin(deg2rad(125.04-1934.136*julianCentury))
%% sunTrueLong - Sun True Long calculated by the calcSunTrueLong function
%% julianCentury - Julian century calculated by the calcJulianCentury function

calcSunAppLong_r(SunTrueLong, JulianCentury) ->
    ?deg2rad(SunTrueLong - 0.00569 - 
		 0.00478*math:sin(?deg2rad(125.04-1934.136*JulianCentury))).

%% Calculate the Mean Obliq Ecliptic in degrees based on the formula: 23+(26+((21.448-julianCentury*(46.815+julianCentury*(0.00059-julianCentury*0.001813))))/60)/60
%% julianCentury - Julian century calculated by the calcJulianCentury function

calcMeanObliqEcliptic(JulianCentury) ->
    (23.0 + (26.0+(21.448-JulianCentury*(46.815+JulianCentury*(0.00059-JulianCentury*0.001813)))/60.0)/60.0).

%% Calculate the Obliq Corr in degrees based on the formula: meanObliqEcliptic+0.00256*cos(deg2rad(125.04-1934.136*julianCentury))
%% meanObliqEcliptic - Mean Obliq Ecliptic calculated by the calcMeanObliqEcliptic function
%% julianCentury - Julian century calculated by the calcJulianCentury function

calcObliqCorr_r(MeanObliqEcliptic, JulianCentury) ->
    ?deg2rad((MeanObliqEcliptic + 
		  0.00256*math:cos(?deg2rad(125.04-1934.136*JulianCentury)))).

%% Calculate the Sun Declination in degrees based on the formula: rad2deg(asin(sin(deg2rad(obliqCorr))*sin(deg2rad(sunAppLong))))
%% obliqCorr - Obliq Corr calculated by the calcObliqCorr function
%% sunAppLong - Sun App Long calculated by the calcSunAppLong function

-spec calcSunDeclination_r(ObliqCorr_r::radians(), SunAppLong_r::radians()) ->
				  radians().

calcSunDeclination_r(ObliqCorr_r, SunAppLong_r) ->
    math:asin(math:sin(ObliqCorr_r) * math:sin(SunAppLong_r)).

%% Calculate the equation of time (minutes) based on the formula:
%% 4*rad2deg(multiFactor*sin(2*deg2rad(geomMeanLongSun))-2*eccentEarthOrbit*sin(deg2rad(geomMeanAnomSun))+4*eccentEarthOrbit*multiFactor*sin(deg2rad(geomMeanAnomSun))*cos(2*deg2rad(geomMeanLongSun))-0.5*multiFactor*multiFactor*sin(4*deg2rad(geomMeanLongSun))-1.25*eccentEarthOrbit*eccentEarthOrbit*sin(2*deg2rad(geomMeanAnomSun)))
%% multiFactor - The Multi Factor vector calculated in the calculate function
%% geomMeanLongSun - The Geom Mean Long Sun vector calculated by the calcGeomMeanLongSun function
%% eccentEarthOrbit - The Eccent Earth vector calculated by the calcEccentEarthOrbit function
%% geomMeanAnomSun - The Geom Mean Anom Sun vector calculated by the calcGeomMeanAnomSun function

-spec calcEquationOfTime(MultiFactor::float(), 
			 GeomMeanLongSun_r::radians(), 
			 EccentEarthOrbit::float(), 
			 GeomMeanAnomSun_r::radians()) -> degree().

%% minutes
calcEquationOfTime(MultiFactor, GeomMeanLongSun_r, 
		   EccentEarthOrbit, GeomMeanAnomSun_r) ->
    Sin_GeomMeanAnomSun_r = math:sin(GeomMeanAnomSun_r),
    A = MultiFactor * math:sin(2.0*GeomMeanLongSun_r),
    B = 2.0 * EccentEarthOrbit * Sin_GeomMeanAnomSun_r,
    C = 4.0 * EccentEarthOrbit * MultiFactor*Sin_GeomMeanAnomSun_r,
    D = math:cos(2.0*GeomMeanLongSun_r),
    E = 0.5 * sqr(MultiFactor)*math:sin(4.0*GeomMeanLongSun_r),
    F = 1.25 * sqr(EccentEarthOrbit)*math:sin(2.0*GeomMeanAnomSun_r),
    4.0 * ?rad2deg(A-B+C*D-E-F).

%% Calculate the HaSunrise in degrees based on the formula: rad2deg(acos(cos(deg2rad(90.833))/(cos(deg2rad(latitude))*cos(deg2rad(sunDeclination)))-tan(deg2rad(latitude))*tan(deg2rad(sunDeclination))))
%% latitude - The latitude defined by the user
%% sunDeclination - The Sun Declination calculated by the calcSunDeclination function

calcHaSunrise(Latitude, SunDeclination_r) ->
    Latitude_r = ?deg2rad(Latitude),
    %% SunDeclination_r = ?deg2rad(SunDeclination),
    ?rad2deg(math:acos(math:cos(?deg2rad(90.833))/(math:cos(Latitude_r)*math:cos(SunDeclination_r)) - math:tan(Latitude_r)*math:tan(SunDeclination_r))).

%% Calculate the Solar Noon based on the formula: (720 - 4 * longitude - equationOfTime + utcOffset * 60) * 60
%% longitude - The longitude is defined by the user
%% equationOfTime - The Equation of Time is calculated by the calcEquationOfTime function
%% utcOffset - The UTC offset is defined by the user

-spec calcSolarNoon(Longitude::degree(), EquationOfTime::degree(),
		    UtcOffset::minute()) -> second().

calcSolarNoon(Longitude, EquationOfTime, UtcOffset) ->
    ((720.0 - 4.0*Longitude - EquationOfTime + UtcOffset*60.0) * 60.0).

%% Compute the number of days between two dates
days_diff(Date1, Date2) ->
    calendar:date_to_gregorian_days(Date2) - calendar:date_to_gregorian_days(Date1).

test() ->
    test(-23.545570, -46.704082, -3.0, {6, 11, 44}, {18, 14, 27}), %% Sao Paulo - Brazil
    test(36.7201600, -4.4203400, 1.0, {7, 16, 45}, {19, 32, 10}),  %% Málaga - Spain
    test(28.613084, 77.209168, 5.5, {6, 21, 45}, {18, 34, 07}),    %% Nova Delhi - India
    test(32.755701, -96.797296, -5.0, {7, 26, 34}, {19, 41, 07}).  %% Dallas - USA

test(Latitude, Longitude, UtcOffset, Sunrise, Sunset) ->
    case get_time(Latitude,Longitude,UtcOffset,{2017,3,23}) of
	{ok,{Sunrise1,Sunset1}} ->
	    Diff1 = abs_time_diff(Sunrise1, Sunrise),
	    Diff2 = abs_time_diff(Sunset1, Sunset),
	    if Diff1 > 0; Diff2 > 0 ->
		    io:format("test fail ~w did not match ~w\n",
			      [{Sunrise1,Sunset1},{Sunrise,Sunset}]),
		    {error, bad_match};
	       true ->
		    ok
	    end;
	{error, Err} ->
	    {error, Err}
    end.

abs_time_diff(Time1, Time0) ->
    abs(time_diff(Time1, Time0)).

time_diff({H,M,S1}, {H,M,S0}) -> (S1-S0);
time_diff({H,M1,S1}, {H,M0,S0}) -> (M1-M0)*60+(S1-S0);
time_diff({H1,M1,S1}, {H0,M0,S0}) -> 
    ((H1-H0)*60*60+(M1-M0)*60+(S1-S0)).

benchmark() ->
    benchmark(100).

benchmark(N) ->
    {Latitude, Longitude, UtcOffset} = eplanet:location(),    
    benchmark(N,Latitude,Longitude,UtcOffset,date()).

benchmark(N,Latitude,Longitude,UtcOffset,Date) ->
    T0 = erlang:monotonic_time(),
    bench_loop(N,Latitude,Longitude,UtcOffset,Date),
    T1 = erlang:monotonic_time(),
    Time = erlang:convert_time_unit(T1 - T0, native, microsecond),
    (N / Time)*1000000.

bench_loop(0,_Latitude,_Longitude,_UtcOffset,_Date) ->
    ok;
bench_loop(I,Latitude,Longitude,UtcOffset,Date) ->
    get_time(Latitude,Longitude,UtcOffset,Date),
    bench_loop(I-1,Latitude,Longitude,UtcOffset,Date).
