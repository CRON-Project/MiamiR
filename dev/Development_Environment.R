
install.packages("renv")
renv::init()        # once
renv::snapshot()   # records the exact environment

git add renv.lock renv
git commit -m "Add renv lockfile"
git push
