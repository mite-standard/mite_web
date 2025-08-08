$(function () {
  $('#builder').queryBuilder({
    plugins: ['bt-tooltip-errors'],
    allow_groups: false,
    allow_empty: true,
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

  $('#builder').queryBuilder('setRules', {
      condition: 'AND',
      rules: []
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
