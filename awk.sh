#!/bin/bash
awk '{t+=$NF} END {print t}'
