%%% @author Tony Rogvall <tony@rogvall.se>
%%% @copyright (C) 2020, Tony Rogvall
%%% @doc
%%%
%%% @end
%%% Created : 12 Apr 2020 by Tony Rogvall <tony@rogvall.se>

-module(xcal).

-compile(export_all).

start() ->
    application:ensure_all_started(erlcron),
    application:ensure_all_started(xbus),
    application:load(xcal).

%% add one job at sunrise and one job at sunset
%% this should be done every day at 00:01

add_sunset_sunrise() ->
    {ok, {Sunrise, Sunset}} = xcal_sun:get_rise_set(),
    io:format("sunrise = ~w\n", [Sunrise]),
    erlcron:cron({{once, Sunrise}, {?MODULE, at_sunrise, []}}),
    io:format("sunset = ~w\n", [Sunset]),
    erlcron:cron({{once, Sunset}, {?MODULE, at_sunset, []}}).


at_sunrise() ->
    xbus:pub("xcal.sunrise", true).

at_sunset() ->
    xbus:pub("xcal.sunset", true).
