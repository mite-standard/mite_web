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
        id: 'orcids',
        label: 'ORCID',
        type: 'string',
        placeholder: 'e.g. 0000-0001-6534-6609',
        validation: {
          format: /^[0-9]/i
        },
        operators: ['equal', 'not_equal', 'contains', 'not_contains']
      },
      {
        id: 'references',
        label: 'Literature DOI',
        type: 'string',
        placeholder: 'e.g. 10.1016/j.chembiol.2020.11.009',
        validation: {
          format: /^10\.\S+$/i
        },
        operators: ['equal', 'not_equal', 'contains', 'not_contains']
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
