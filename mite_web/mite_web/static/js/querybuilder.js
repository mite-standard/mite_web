$(function () {
  $('#builder').queryBuilder({
    plugins: ['bt-tooltip-errors'],
    filters: [
      {
        id: 'accession',
        label: 'Accession',
        type: 'string'
      },
      {
        id: 'contributor',
        label: 'Contributor',
        type: 'string'
      },
      {
        id: 'reviewer',
        label: 'Reviewer',
        type: 'string'
      },
    ]
  });

  $('#btn-get').on('click', function () {
    const result = $('#builder').queryBuilder('getRules');
    if (!result) {
      $('#output').text('Invalid query');
    } else {
      $('#output').text(JSON.stringify(result, null, 2));
    }
  });
});
