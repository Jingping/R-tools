#!/bin/bash

git init
touch README
git add README
git commit -m 'first commit'
git remote add origin git@github.com:Jingping/R-tools.git
git push -u origin master
