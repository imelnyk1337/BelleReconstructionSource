scp -p reco.* userinfo.* imelnyk@sshcc1.kek.jp:
ssh -t imelnyk@sshcc1.kek.jp scp -p reco.* userinfo.* login.cc.kek.jp:DSDS2317/reco/src
ssh -t imelnyk@sshcc1.kek.jp ssh login.cc.kek.jp ./cmd_make.sh
