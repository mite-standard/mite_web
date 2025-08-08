$(function () {
  $('#builder').queryBuilder({
    plugins: ['bt-tooltip-errors'],
    allow_groups: false,
    allow_empty: true,
    operators: [
      'equal',
      'not_equal',
      'contains',
      'not_contains',
      'is_null',
      'is_not_null'
    ],
    filters: [
      {
        id: 'accession',
        label: 'Accession',
        type: 'string'
      },
      {
        id: 'persons.orcid',
        label: 'ORCID',
        type: 'string'
      },




      {
        id: 'enzyme.mibig_id',
        label: 'MIBiG ID',
        type: 'string'
      },
      {
        id: 'enzyme.uniprot_id',
        label: 'UniProt ID',
        type: 'string'
      },
      {
        id: 'enzyme.cofactors.cofactor_name',
        label: 'Cofactor name',
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
