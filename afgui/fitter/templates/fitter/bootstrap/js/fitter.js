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
    var retval = {};
    var count = parseInt($('#num_dragpoints').val());
    for (var i=1; i<=count; i++) {
        x = parseInt($('#' + which + '_dragline' + i).css('left').slice(0,-2)) + 5;
        retval['x'+i] = x;
        y = parseInt($('#' + which + '_dragline' + i).css('top').slice(0,-2)) + 5;
        retval['y'+i] = y;
    }
    return retval
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

function update_draglines() {
    coords = get_clean_coords('top');
    var line = '';
    var count = parseInt($('#num_dragpoints').val());
    for (var i=1; i<=count; i++) {
        line += (i==1 ? 'M ' : 'L ') + coords['x'+i] + ' ' + coords['y'+i];
    }
    $('#top_dragline').attr('d',line);
    coords = get_clean_coords('bottom');
    var line = '';
    for (var i=1; i<=count; i++) {
        line += (i==1 ? 'M ' : 'L ') + coords['x'+i] + ' ' + coords['y'+i];
    }
    $('#bottom_dragline').attr('d',line);
}

function create_draglines(count) {
    var botdragline_top = (height - 5) + 'px';
    $('#top_dragpoints').html('');
    $('#bottom_dragpoints').html('');
    for (var i=1; i<=count; i++) {
        $('<div/>', {
            id: 'top_dragline' + i,
            class: "smallcirc",
        }).appendTo($('#top_dragpoints')).css({"background-color": "darkgreen",
                                               "left":parseInt(width*(i-1)/(count-1))+"px",
                                               "top":"0px"});
        $('<div/>', {
            id: 'bottom_dragline' + i,
            class: "smallcirc",
        }).appendTo($('#bottom_dragpoints')).css({"background-color": "darkred",
                                               "left":parseInt(width*(i-1)/(count-1))+"px",
                                               "top":botdragline_top});
    }
    // update coordinates
    update_draglines();
    // make draggable
    $('.smallcirc').draggable({
        containment: "#fig_h_orig",
        drag: function() {
            update_draglines();
        },
        stop: function() {
            var max_left = parseInt($('#fig1').attr('data-mywidth'));
            var left = $(this).css('left');
            if (left.slice(0,-2) > max_left) {
                $(this).css('left',max_left+'px');
                update_draglines();
            }
        }
    });
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
        var count = parseInt($('#num_dragpoints').val());
        create_draglines(count);
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

function get_formatted_coords(which) {
    coords = get_clean_coords(which);
    var text_coords = '';
    var count = parseInt($('#num_dragpoints').val());
    for (var i=1; i<=count; i++) {
        text_coords += '(' + coords['x'+i] + ',' + coords['y'+i] + '),';
    }
    return text_coords.slice(0,-1);
}

function clean_data(callback, threshold) {
    var text_top_coords = get_formatted_coords('top');
    var text_bottom_coords = get_formatted_coords('bottom');
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

