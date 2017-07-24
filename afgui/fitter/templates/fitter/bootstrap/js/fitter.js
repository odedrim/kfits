function input_cleanup() {
    $('#input_path').css('background-color','#FFFFFF');
    clear_user_message();
}

function fig_cleanup() {
    $('#fig1').html('');
    $('#fig2').html('');
    $('#fig3').html('');
    $('#fig4').html('');
    $('.smallcirc').hide();
    $('#top_dragline').attr('d','');
    $('#bottom_dragline').attr('d','');
    $('#approx_start').hide();
    $('#fit_data').hide();
    $('#gross_filter').hide();
    $('.auto_kinetics').remove();
    $('#model_name').html('');
    $('#t1').html('');
    $('#t2').html('');
    set_progress($('#total_progress'), 0, '');
}

function get_clean_coords(which) {
    x1 = parseInt($('#' + which + '_dragline1').css('left').slice(0,-2)) + 5;
    y1 = parseInt($('#' + which + '_dragline1').css('top').slice(0,-2)) + 5;
    x2 = parseInt($('#' + which + '_dragline2').css('left').slice(0,-2)) + 5;
    y2 = parseInt($('#' + which + '_dragline2').css('top').slice(0,-2)) + 5;
    x3 = parseInt($('#' + which + '_dragline3').css('left').slice(0,-2)) + 5;
    y3 = parseInt($('#' + which + '_dragline3').css('top').slice(0,-2)) + 5;
    return {x1:x1, y1:y1, x2:x2, y2:y2, x3:x3, y3:y3};
}

function set_progress(bar, value, text) {
    bar.css('width', value+'%');
    bar.attr('aria-valuenow', value);
    bar.html(text);
}

function check_input_file(fname, callback) {
    var res = $.get("backend", {"function":"check_input_file", "fnames":fname}, function(data) {
        callback(data[0], data[1]);
    });
}

function is_ok_to_run() {
    // return $('#input_path').css('background-color') == "rgb(238, 250, 245)";
    return $('#fdata').hasClass("btn-success");
}

function message_user(content, alert_type) {
    if (alert_type == null) {
        alert_type = "success";
    }
    $('#input_message_user').html(content);
    $('#input_message_user').attr('class', 'alert alert-'+alert_type);
}

function clear_user_message() {
    $('#input_message_user').html("");
    $('#input_message_user').attr('class', 'hidden');
}

function run_analysis() {
    // cleanup everything
    fig_cleanup();
    // load figure
    var img = $('<img />', {attr: {src: 'backend?function=plot_data&fnames=' + $('#input_path').val() + '&_=' + new Date().getTime()}}).load(function() { 
        // get figure size
        width = parseInt($('#fig1').attr('data-mywidth'));
        height = parseInt($('#fig1').attr('data-myheight'));
        // create draglines
        $('#top_dragline1').css({"left":"0px", "top":"0px"}).show();
        $('#top_dragline2').css({"left":parseInt(width/2)+"px", "top":"0px"}).show();
        $('#top_dragline3').css({"left":width+"px", "top":"0px"}).show();
        coords = get_clean_coords('top');
        $('#top_dragline').attr('d','M ' + coords.x1 + ' ' + coords.y1 + ' L ' + coords.x2 + ' ' + coords.y2 + ' L ' + coords.x3 + ' ' + coords.y3);
        var botdragline_top = (height - 5) + 'px';
        $('#bottom_dragline1').css({"left":"0px", "top":botdragline_top}).show();
        $('#bottom_dragline2').css({"left":parseInt(width/2)+"px", "top":botdragline_top}).show();
        $('#bottom_dragline3').css({"left":width+"px", "top":botdragline_top}).show();
        coords = get_clean_coords('bottom');
        $('#bottom_dragline').attr('d','M ' + coords.x1 + ' ' + coords.y1 + ' L ' + coords.x2 + ' ' + coords.y2 + ' L ' + coords.x3 + ' ' + coords.y3);
        // move to original curve view
        $('#fig_b_orig').click();
        $('#gross_filter').show();
        // progress
        set_progress($('#total_progress'), 33.33, 'Data Loaded');
    })
    .each(function() {
        //Cache fix for browsers that don't trigger .load() - https://stackoverflow.com/questions/2392410/jquery-loading-images-with-complete-callback
        if(this.complete) $(this).trigger('load');
    });
    img.appendTo($('#fig1'));
}

function clean_data(callback, threshold) {
    coords = get_clean_coords('top');
    var text_top_coords = '(' + coords.x1 + ',' + coords.y1 + '),(' + coords.x2 + ',' + coords.y2 + '),(' + coords.x3 + ',' + coords.y3 + ')';
    coords = get_clean_coords('bottom');
    var text_bottom_coords = '(' + coords.x1 + ',' + coords.y1 + '),(' + coords.x2 + ',' + coords.y2 + '),(' + coords.x3 + ',' + coords.y3 + ')';
    if (typeof threshold == 'undefined') {
        var func = 'clean_data_optimise_noise_threshold';
    } else {
        var func = 'clean_data&noise_threshold=' + threshold;
    }
    $.getJSON('backend?function=' + func + '&fnames=' + $('#input_path').val() + '&model=' + $('#model_choice').val() + '&threshold_points=' + text_top_coords + '&rev_threshold_points=' + text_bottom_coords + '&approx_start=' + $('#approx_start').css('left').slice(0,-2) + ($('#noise_only_above').is(':checked') ? '&noise_only_above' : '') + ($('#reaches_plateau').is(':checked') ? '&search_for_end' : ''), function(data) {
        // cleanup
        $('.auto_kinetics').remove();
        // parse data
        d = data[0];
        $('#fig4').html(d[0]);
        $('#threshold').val(d[1]);
        // print kinetic parameters
        $('#model_name').html(d[2]);
        $('#t1').html(d[3]);
        $('#t2').html(d[4]);
        for(var key in d[5]){
            if (d[5].hasOwnProperty(key)) {
                var tr = $('<tr></tr>').addClass('auto_kinetics');
                tr.append($('<td></td>').html(key));
                tr.append($('<td></td>').html(d[5][key]));
                $('#kinetic_variables').append(tr);
            }
        }
        callback(text_top_coords, text_bottom_coords);
    });
}

