#!/bin/zsh
# Double-click to start the MsiaGen UI. Keep this window open while using it.
cd "$(dirname "$0")"
clear
echo ""
echo "  =========================================================== "
echo "  MSIAGEN is starting. Please wait, this takes about 1 minute."
echo "  Your web browser will open by itself when it is ready."
echo "  Do not close this window while you use MsiaGen."
echo "  =========================================================== "
echo ""
julia --project=. -t auto start.jl
