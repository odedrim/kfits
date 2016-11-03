function check_input_file(fname, callback) {
    var res = $.get("backend", {"function":"check_input_file", "fnames":fname}, function(data) {
        callback(data[0], data[1]);
    });
}
