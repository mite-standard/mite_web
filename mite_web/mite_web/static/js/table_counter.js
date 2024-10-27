window.onload = function () {
    $('#table').on('search.bs.table', function () {
        const searchResults = $('#table').bootstrapTable('getData');
        console.log(searchResults.length);
        updateSearchCount(searchResults.length);
    });

    function updateSearchCount(count) {
        $('#search-count').text(`Rows found: ${count}`);
    }
};