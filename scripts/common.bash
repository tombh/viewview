# shellcheck disable=2120
function _panic() {
	local message=$1
	echo >&2 "$message"
	exit 1
}

function _pushd {
	# shellcheck disable=2119
	command pushd "$@" >/dev/null || _panic
}

function _popd {
	# shellcheck disable=2119
	command popd "$@" >/dev/null || _panic
}
