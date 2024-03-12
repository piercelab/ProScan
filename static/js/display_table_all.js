$(document).ready(function(){
            $.ajax({
              url: '/fetch_table/'+ encodeURIComponent(jobName),
              type: 'GET',
              success: function(response){
                $('#table_display').append(response.table);
                $('#table_display table').tablesorter({
                theme: 'default',
                sortReset: true,
                widgets: ['filter','stickyHeaders'], 
                widgetOptions: {                
      // jQuery selector or object to attach sticky header to
      stickyHeaders_attachTo : '#table_display', 
      scrollableArea: '#table_display'}
                });     
              },
                error: function(response){
                  console.log('Error fetching table');
                }
            });
          });
