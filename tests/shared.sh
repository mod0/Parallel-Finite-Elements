projectDir="../"
executable="ellipticmain"

function timeRun {
    (time ${projectDir}${executable} $1 $2 $3 &>/dev/null) 2>&1 | grep real
}

function build {
    make f=$1 -C ${projectDir} &>/dev/null
}
