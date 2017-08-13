$(document).ready(function() {
    var MAX_UPLOAD_SIZE = 10*1024*1024;
    $.ajaxSetup({ cache: false });

    $.get('backend?function=get_email_oded', function(data) {
        $('#email_oded').html(data.replace('123','@'));
    })
    $.get('backend?function=get_email_dana', function(data) {
        $('#email_dana').html(data.replace('123','@'));
    })

    $(document).keyup(function(e){
        // escape
        if(e.keyCode === 27) {
            $('#m_home').click();
        }
    });

    $('#num_dragpoints').change(function() {
        if ($(this).val() > 9) {
            $(this).val(9);
        } else if ($(this).val() < 2) {
            $(this).val(2);
        }
        // recreate draglines iff they are already visible (otherwise they will be created later, no worries...)
        if (($('#top_dragline1').length > 0) && ($('#top_dragline1').css('display') != "none")) {
            create_draglines($(this).val());
        }
    });

    $('#approx_start').draggable({
        stop: function() {
            var max_left = parseInt($('#fig1').attr('data-mywidth'));
            $(this).css('top','0px');
            var left = $(this).css('left');
            if (left.substr(0,1) == '-') {
                $(this).css('left','0px');
            } else if (left.slice(0,-2) > max_left) {
                $(this).css('left',max_left+'px');
            }
        }
    });

    $(':file').change(function() {
        var file = this.files[0];
        var size = file.size;
        var type = file.type;
        if (size > MAX_UPLOAD_SIZE) {
            message_user("<strong>Notice!</strong> You are about to upload a file bigger than what the system can handle!", "warning");
            $("#fdata").attr("class", "btn btn-warning");
        } else if ((type != 'text/plain') && (type != 'text/csv') && (type != 'application/vnd.ms-excel')) {
            if (type == "") {
                type = "UNKNOWN FILETYPE";
            }
            message_user("<strong>Notice!</strong> The file you are about to upload is in a format the system cannot handle ("+type+")", "warning");
            $("#fdata").attr("class", "btn btn-warning");
        } else {
            clear_user_message();
            $("#fdata").attr("class", "btn btn-primary");
        }
    });

    $('#upload').click(function() {
            var formData = new FormData($('form')[0]);
            $.ajax({
                url: 'upload_text.htm',
                type: 'POST',
                xhr: function() {  // Custom XMLHttpRequest
                    var myXhr = $.ajaxSettings.xhr();
                    if (myXhr.upload) { // Check if upload property exists
                        myXhr.upload.addEventListener('progress', function(e) {
                            if(e.lengthComputable){
                                //$('').attr({value:e.loaded,max:e.total});
                                set_progress($('#total_progress'), 33.33*e.loaded/e.total, 'Data Loading');
                            }
                        }, false); // For handling the progress of the upload
                    }
                    return myXhr;
                },
                //Ajax events
                beforeSend: function() {},
                success: function(data) {
                    if (data != false) {
                        $("#input_path").val(data);
                        check_input_file(data, function(is_good, reason) {
                            if (is_good) {
                                $('#fdata').attr('class', 'btn btn-success');
                                run_analysis();
                            } else {
                                $('#fdata').attr('class', 'btn btn-danger');
                                message_user("<strong>Cannot Load File!</strong> " + reason, "danger");
                            }
                        });
                    }
                },
                error: function() { console.log("ERROR"); },
                // Form data
                data: formData,
                //Options to tell jQuery not to process data or worry about content-type.
                cache: false,
                contentType: false,
                processData: false
            });
    });

    $('.loadExample').click(function() {
        $('#fdata').attr('class', 'btn btn-success');
        $("#input_path").val($(this).attr('path'));
        run_analysis();
    });

    $('#input_clear').bind('click', function() {
        $('#fdata').val("");
        $('#fdata').attr("class", "btn btn-primary");
        $('#input_path').val("");
        input_cleanup();
        fig_cleanup();
    });

    $('#input_path').on('change keypress keyup paste', function() {
        input_cleanup();
    });

    $('#model_choice').change(function() {$('#gross_filter').trigger('click')});
    $('#noise_only_above').change(function() {$('#gross_filter').trigger('click')});
    $('#reaches_plateau').change(function() {$('#gross_filter').trigger('click')});

    $('#gross_filter').bind('click', function() {
        if (is_ok_to_run()) {
            var text_top_coords = get_formatted_coords('top');
            var text_bottom_coords = get_formatted_coords('bottom');
            $('#fig2').html('<img src="backend?function=plot_data&fnames=' + $('#input_path').val() + '&threshold_points=' + text_top_coords + '&rev_threshold_points=' + text_bottom_coords + '&_=' + new Date().getTime() + '">');
            $('#approx_start').show();
            $('#fig_b_filt').click();
            $('#fit_data').show();
            set_progress($('#total_progress'), 66.67, 'Outliers Filtered');
        }
    });

    $('.fig_badge').bind('click', function() {
        var this_id = $(this)[0].id;
        var head_id = '#fig_h_' + this_id.substr(6);
        $('.fig_badge').removeClass('active');
        $('.fig_head').hide();
        $(this).addClass('active');
        $(head_id).show();
    });

    $('#fit_data').bind('click', function() {
        // clean data with automatic noise threshold optimisation
        clean_data(function(text_top_coords, text_bottom_coords) { 
            $('#fig3').html('<img src="backend?function=fit_data&fnames=' + $('#input_path').val() + '&model=' + $('#model_choice').val() + '&threshold_points=' + text_top_coords + '&rev_threshold_points=' + text_bottom_coords + '&approx_start=' + $('#approx_start').css('left').slice(0,-2) + ($('#reaches_plateau').is(':checked') ? '&search_for_end' : '') + '">');
            $('#fig_b_fit').click();
            $('#arrow_to_results').delay(2000).fadeIn(1000).fadeOut(1000).fadeIn(1000).fadeOut(1000).fadeIn(1000).fadeOut(1000);
            set_progress($('#total_progress'), 100, 'Fitting Done!');
        });
    });

    $('.topnav_btn').bind('click', function() {
        $('.topnav_btn').removeClass('active');
        $(this).addClass('active');
        if ($(this)[0] == $('#m_home')[0]) {
            $('.overall').hide();
        } else {
            $('.overall').show();
        }
        $('.message').hide();
        if ($(this)[0] == $('#m_instructions')[0]) {
            $('#instructions').show();
        } else if ($(this)[0] == $('#m_about')[0]) {
            $('#about').show();
        } else if ($(this)[0] == $('#m_contact')[0]) {
            $('#contact').show();
        }
    });

    $('.navbar-brand').bind('click', function() {
        $('#m_home').click();
    });
    $('.dropdown').bind('click', function() {
        $('.overall').hide();
        $('.topnav_btn').removeClass('active');
        $('#m_home').addClass('active');
    });

    $('#threshold').on('input', function() {
        $('#clean_data').attr('disabled', false);
    });

    $('#clean_data').bind('click', function() {
        if (!$(this).hasClass('grayed_out_blue')) {
            clean_data(function() {
                $('#clean_data').attr('disabled', true);
            }, $('#threshold').val());
        }
    });

    $('#get_clean_data').bind('click', function() {
        var text_top_coords = get_formatted_coords('top');
        var text_bottom_coords = get_formatted_coords('bottom');
        if ($('#fdata').get(0).files.length) {
            var filename = $('#fdata').get(0).files[0].name;
        } else {
            var filename = 'unknown';
        }
        var sid = $('#input_path').val().slice(-36,-4);
        $('#dlink').html('Session ID: <a href=get.htm?sid=' + sid + '>' + sid + '</a> (Valid for 24 hours)');
        window.location = 'backend?function=get_clean_data&fnames=' + $('#input_path').val() + '&model=' + $('#model_choice').val() + '&noise_threshold=' + $('#threshold').val() + '&threshold_points=' + text_top_coords + '&rev_threshold_points=' + text_bottom_coords + '&approx_start=' + $('#approx_start').css('left').slice(0,-2) + '&output_fnames=' + filename + '_clean.csv' + ($('#noise_only_above').is(':checked') ? '&noise_only_above' : '') + ($('#reaches_plateau').is(':checked') ? '&search_for_end' : '');
    });
});
