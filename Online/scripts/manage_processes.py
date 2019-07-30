#!/usr/bin/env python3

"""
Control script to manage background processes such as
runinfod and related epics IOC, and auto-replay.

Author: Sylvester Joosten <sjjooste@jlab.org>

Commands:
    start
    stop
    restart
    expert
"""

import os
import libtmux
import argparse
import time

TMUX_SOCKET_NAME = 'vcs_online'
TMUX_SESSION_NAME = 'background'
EPICS_DIR = '/u/group/c-csv/whit/software/epics/apps/hallc_epics_run_info/iocBoot/ioctest1'
EPICS_CMD = './st.cmd'
RUNINFOD_CMD = 'bash Online/scripts/launch_runinfod'
AUTOREPLAY_CMD = 'python Online/scripts/start_auto_replay.py'

class ProcessAlreadyRunningError(Exception):
    pass
class NoProcessRunningError(Exception):
    pass
class UnknownArgumentError(Exception):
    pass

def _replay_dir():
    """Get the current replay directory.

    Current replay directory is the replay directory that contains this script.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.realpath(script_dir + "/../..")

def start():
    """Start the background processes.

    Will raise an exception if already running.
    """
    server = libtmux.Server(socket_name=TMUX_SOCKET_NAME)
    session = server.find_where({"session_name": TMUX_SESSION_NAME})
    ## Ensure the session is not currently running
    if session:
        raise ProcessAlreadyRunningError()
    print("Starting background session...")
    session = server.new_session(
            session_name=TMUX_SESSION_NAME, 
            start_directory=_replay_dir())
    ## Start the epics IOC
    print("Starting EPICS IOC...")
    epics_win = session.new_window(
            window_name="epics",
            start_directory=EPICS_DIR,
            attach=False,
            window_shell=EPICS_CMD)
    print("   ...Waiting 10s")
    time.sleep(10)
    ## Start runinfod
    print("Starting runinfod...")
    runinfod_win = session.new_window(
            window_name="runinfod",
            start_directory=_replay_dir(),
            attach=False,
            window_shell=RUNINFOD_CMD)
    print("Background processess started.")
    ## Start auto replay
    print("Starting auto-replay...")
    replay_win = session.new_window(
            window_name="replay_manager",
            start_directory=_replay_dir(),
            attach=False,
            window_shell=AUTOREPLAY_CMD)

def stop(restart=False):
    """Stop the background processes.

    Will raise an exception if no processes are running.

    Argument:
        * restart: if set to True we are in restart mode, no
                   exception will be raised if the the process is not
                   running.
    """
    print("Stopping background processes...")
    server = libtmux.Server(socket_name=TMUX_SOCKET_NAME)
    session = server.find_where({"session_name": TMUX_SESSION_NAME})
    ## Ensure the session is running, except when in restart mode
    if not session:
        if restart:
            print("...Nothing to be done.")
            return
        raise NoProcessRunningError()
    print("Killing the background session...")
    session.kill_session()
    print("Background processes stopped.")

def restart():
    """Restart the background processes.

    Also works when no processes are running.
    """
    print("Restarting background processes.")
    stop(restart=True)
    time.sleep(1)
    start()

def expert():
    """Expert mode: get manual control on tmux session."""
    server = libtmux.Server(socket_name=TMUX_SOCKET_NAME)
    session = server.find_where({"session_name": TMUX_SESSION_NAME})
    print("Starting expert mode... attaching to session", TMUX_SESSION_NAME)
    ## Ensure the session is running
    if not session:
        raise NoProcessRunningError()
    session.attach_session()
    print("Expert mode done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "command",
            help="Known commands are start, stop, restart and expert")
    args = parser.parse_args()
    if args.command not in ['start', 'stop', 'restart', 'expert']:
        parser.print_help()
        raise UnknownArgumentError
    locals()[args.command]()
