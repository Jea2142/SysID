#!/bin/sh

# Set Session Name
SESSION="FlightRecord"
SESSIONEXISTS=$(tmux list-sessions | grep $SESSION)

# Only create tmux session if it doesn't already exist
if [ "$SESSIONEXISTS" = "" ]
then
	# Start New Session with our name
	tmux new-session -d -s $SESSION

	# Name first Pane and split into multiple
	tmux rename-window -t 0 'Gazebo Tmux Test'

	tmux select-pane -t 0
	tmux split-window -h   




	# Send commands to pane 0
	NAME="Zara Ali"
	echo $NAME
	tmux send-keys -t 0 'sleep 5' C-m
	tmux send-keys -t 0 'cd ~/catkin_ws/src' C-m 


	# Send commands to pane 1
	tmux send-keys -t 1 'sleep 5' C-m
	tmux send-keys -t 1 'cd ~/src/Firmware' C-m






fi

# Attach Session, on the Main window
tmux attach-session -t $SESSION:0
