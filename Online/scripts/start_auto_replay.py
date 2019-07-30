#!/usr/bin/env python3

"""
Control to launch auto-replay jobs for VCS as a background process.

Note that this script is meant to be controlled indirectly through
manage_processes.py, with which it shares a tmux socket

Author: Sylvester Joosten <sjjooste@jlab.org>

"""

from hallc.run_daemon import RunDaemon
import os
import sys
import libtmux
import json

TMUX_SOCKET = 'vcs_online'
TMUX_SESSION = 'replay'
REPLAY_CMD = './run_full_replay.csh {} {} {}'
DBDIR='Online/database'

def _replay_dir():
    """Get the current replay directory.

    Current replay directory is the replay directory that contains this script.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.realpath(script_dir + "/../..")

def _get_analysis_command(run):
    """Get the analysis command for this run"""
    run = str(run)
    print('Launching full replay for run', run)
    runinfo = {}
    print('  --> Loading run info')
    with open('{}/rundb_coin.json'.format(DBDIR)) as f:
        rundb = json.load(f)
        runinfo = rundb[run]
    target = runinfo['target']['target_label']
    print('  --> Target set to:', target)
    print('  --> Loading spectrometer info...')
    hms_th = runinfo['spectrometers']['hms_angle']
    hms_p = runinfo['spectrometers']['hms_momentum']
    shms_th = runinfo['spectrometers']['shms_angle']
    shms_p = runinfo['spectrometers']['shms_momentum']
    print('      + HMS: {} GeV and {} degrees'.format(hms_p, hms_th))
    print('      + SHMS: {} GeV and {} degrees'.format(hms_p, hms_th))
    kin = 'other'
    with open('{}/settings.json'.format(DBDIR)) as f:
        settings = json.load(f)
        for name in settings:
            setting = settings[name]
            if setting['target'] == target and \
                    abs(setting['hms_angle'] - hms_th) < 0.1 and \
                    abs(setting['hms_momentum'] - hms_p) < 0.015 and \
                    abs(setting['shms_angle'] - shms_th) < 0.1 and \
                    abs(setting['shms_momentum'] - shms_p) < 0.015:
                kin = setting['name']
    print('  --> Kinematic identified as:', kin)
    ## run either VCS or elastics
    ana_opt = 'vcs' if 'Kin' in kin else 'ela'
    ## dummy has to be lower case for analysis script
    if (target == 'DUMMY'):
        target = target.lower()
    print('  --> Analysis options:', run, ana_opt, target)
    return REPLAY_CMD.format(run, ana_opt, target)

def launch(run):
    """Launch replay for this run as a tmux session."""
    ## Get our replay command
    command = _get_analysis_command(run)
    print('Replay command: "{}"'.format(command))

    print('Creating new tmux window for the replay')
    server = libtmux.Server(socket_name=TMUX_SOCKET)
    session = server.find_where({"session_name": TMUX_SESSION})
    ## Ensure the session is running
    if not session:
        session = server.new_session(
                session_name=TMUX_SESSION,
                start_directory=_replay_dir())
    ## create a new window for this run
    window_name = "run_{}".format(run)
    window = session.find_where({"window_name": window_name})
    if window:
        window.kill_window()
    print('Launching replay')
    window = session.new_window(
            window_name=window_name,
            start_directory=_replay_dir(),
            attach=False,
            window_shell=command)

if __name__ == "__main__":
    ## Attach the launch(run) action to the run_stop event for the COIN DAQ
    d = RunDaemon()
    d.on_event('run_stop', launch, run_type='coin')
    d.start()
