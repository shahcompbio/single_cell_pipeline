CURR_HEAD=$(git rev-parse $(git rev-parse --abbrev-ref HEAD))
TAG=$(git describe --tags $(git rev-list --tags --max-count=1))
TAG_HEAD=$(git rev-parse $TAG^{commit})

if test $CURR_HEAD != $TAG_HEAD; then
    echo "Branch is not tagged"
    exit -1
fi

