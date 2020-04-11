%%% @author Tony Rogvall <tony@rogvall.se>
%%% @copyright (C) 2020, Tony Rogvall
%%% @doc
%%%
%%% @end
%%% Created : 12 Apr 2020 by Tony Rogvall <tony@rogvall.se>

-module(sunset_cron).

-compile(export_all).

start() ->
    erlcron:start(),
    xbus:pub_meta("xcal.sunrise", [{unit,"bool"}]),
    xbus:pub_meta("xcal.sunset", [{unit,"bool"}]),
    erlcron:cron(setup_sunset_sunrise,
		 {{daily, {00,01,00}}, 
		  {?MODULE, add_sunset_sunrise, []}}).

%% add one job at sunrise and one job at sunset
%% this should be done every day at 00:01

add_sunset_sunrise() ->
    {ok, Sunset, Sunrise} = sunset_sunrise:get_time(),
    erlcron:cron(sunrise, {once, Sunrise, {?MODULE, at_sunrise, []}}),
    erlcron:cron(sunset, {once, Sunset, {?MODULE, at_sunset, []}}).

at_sunrise() ->
    xbus:pub("xcal.sunrise", true).

at_sunset() ->
    xbus:pub("xcal.sunset", true).
