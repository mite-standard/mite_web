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

  $('form').on('submit', function (e) {
    const rules = $('#builder').queryBuilder('getRules');
    if (!rules || !rules.valid) {
      alert("Invalid query!");
      e.preventDefault(); // stop form from submitting
      return;
    }
    $('#rules-input').val(JSON.stringify(rules));
  });
});
