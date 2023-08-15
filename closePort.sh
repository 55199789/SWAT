sudo kill -9 $(pidof bin/serverProxy)
sudo kill -9 $(pidof bin/clientProxy)
sudo kill -9 $(lsof -t -i :12341)
# ./build.sh