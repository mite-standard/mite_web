// Checks for Firefox Mobile browser which is not supported

function testFirefoxMobile() {
    if (navigator.userAgent.indexOf("Firefox") != -1 && /Mobile/i.test(navigator.userAgent)) {
        alert("Sorry, form completion with Firefox Mobile is not supported on this site. Please use an alternative browser or use Firefox Desktop");
    };
}

// prevents repeated POST requests
function disableButton(form, button) {
  form.addEventListener('submit', function (event) {
      if (!form.checkValidity()) {
          return;
      }
      setTimeout(() => {
            button.disabled = true;
            button.innerText = 'Submitting...';
        }, 0);
  });
}


// Removes a Div of a specified class closest to the button; used by all elements that can be added/removed by user-input
function removeField(button, class_ref) {
    const targetDiv = button.closest(class_ref);
    if (targetDiv) {
        targetDiv.remove();
    }
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeName(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_name" name="enzyme_name" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.enzyme?.name ?? ""}' required title="The commonly used enzyme name (e.g. 'McjB')">
            <label for="enzyme_name" class="form-label" >Name<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeDescription(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_description" name="enzyme_description" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.enzyme?.description ?? ""}' maxlength="50" required title="Enzyme type description (e.g. 'O-methyltransferase')">
            <label for="enzyme_description" class="form-label">Type<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createGeneralComment(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="comment" maxlength="300" name="comment" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.comment ?? ""}' title="Brief free-text comment on enzyme and/or reaction details not applicable to other form fields.">
            <label for="comment" class="form-label">General comment <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeUniprot(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_uniprot" name="enzyme_uniprot" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.enzyme?.databaseIds?.uniprot ?? ""}' title="UniProtKB or UniParc ID of enzyme">
            <label for="enzyme_description" class="form-label" >UniProt ID<span class="mandatory">**</span> <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeGenpept(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_genpept" name="enzyme_genpept" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.enzyme?.databaseIds?.genpept ?? ""}' title="NCBI GenPept ID of enzyme">
            <label for="enzyme_genpept" class="form-label">GenPept ID<span class="mandatory">**</span> <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeMibig(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_mibig" name="enzyme_mibig" pattern="^BGC.+" minlength="10" maxlength="10" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.enzyme?.databaseIds?.mibig ?? ""}' title="MIBiG ID of the BGC containing the enzyme">
            <label for="enzyme_mibig" class="form-label">MIBiG ID <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeWikidata(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_wikidata" name="enzyme_wikidata" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.enzyme?.databaseIds?.wikidata ?? ""}' title="The Wikidata QID of the enzyme">
            <label for="enzyme_wikidata" class="form-label">Wikidata ID <i class="bi bi-question-circle"></i></label>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeCofactors(container, data = {}, form_vals) {
    const entryHtml = `
        <div class="card card-body">
            <div class="row d-flex align-items-center">
                <div class="col">
                    <div class="form-label mb-2">Inorganic cofactors <i class="bi bi-question-circle" data-bs-toggle="tooltip" title="Inorganic cofactors required for the functioning of the enzyme"></i></div>
                </div>
                <div class="col-auto">
                    <button class="accordion-button-info collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#enzyme-cofactors-inorganic" aria-controls="enzyme-cofactors-inorganic" aria-expanded="false"></button>
                </div>
            </div>
            <div class="row collapse" id="enzyme-cofactors-inorganic"></div>
            <div class="row d-flex align-items-center mt-2">
                <div class="col">
                    <div class="form-label mb-2">Organic cofactors <i class="bi bi-question-circle" data-bs-toggle="tooltip" title="Organic cofactors required for the functioning of the enzyme"></i></div>
                </div>
                <div class="col-auto">
                    <button class="accordion-button-info collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#enzyme-cofactors-organic" aria-controls="enzyme-cofactors-organic" aria-expanded="false"></button>
                </div>
            </div>
            <div class="row collapse" id="enzyme-cofactors-organic"></div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);

    const containerInorganic = document.getElementById(`enzyme-cofactors-inorganic`);
    for (let inorg of form_vals.inorganic) {
        if (data.enzyme.cofactors && data.enzyme.cofactors.inorganic && data.enzyme.cofactors.inorganic.includes(inorg) ) {
            const checkHtmlInorg = document.createElement('div');
            checkHtmlInorg.classList.add('col-1');
            checkHtmlInorg.classList.add('form-text');
            checkHtmlInorg.innerHTML = addTickedCheckbox(`enzyme-cofactors-inorganic-check[]`, inorg);
            containerInorganic.appendChild(checkHtmlInorg)
        } else {
            const checkHtmlInorg = document.createElement('div');
            checkHtmlInorg.classList.add('col-1');
            checkHtmlInorg.classList.add('form-text');
            checkHtmlInorg.innerHTML =  addUntickedCheckbox(`enzyme-cofactors-inorganic-check[]`, inorg);
            containerInorganic.appendChild(checkHtmlInorg)
        }
    }

    const containerOrganic = document.getElementById(`enzyme-cofactors-organic`);
    for (let org of form_vals.organic) {
        if (data.enzyme.cofactors && data.enzyme.cofactors.organic && data.enzyme.cofactors.organic.includes(org) ) {
            const checkHtmlOrg = document.createElement('div');
            checkHtmlOrg.classList.add('col-3');
            checkHtmlOrg.classList.add('text-break');
            checkHtmlOrg.classList.add('form-text');
            checkHtmlOrg.innerHTML = addTickedCheckbox(`enzyme-cofactors-organic-check[]`, org);
            containerOrganic.appendChild(checkHtmlOrg)
        } else {
            const checkHtmlOrg = document.createElement('div');
            checkHtmlOrg.classList.add('col-3');
            checkHtmlOrg.classList.add('text-break');
            checkHtmlOrg.classList.add('form-text');
            checkHtmlOrg.innerHTML =  addUntickedCheckbox(`enzyme-cofactors-organic-check[]`, org);
            containerOrganic.appendChild(checkHtmlOrg)
        }
    }
}


// Populate the enzyme reference forms, triggered by DOMContentLoaded event
function populateEnzymeRefForm(container, data) {
    if (data.enzyme?.references) {
        data.enzyme.references.forEach(entry => {
            const entryHtml = createHtmlRef(entry, "enzyme_ref", "enzyme_ref[]");
            container.insertAdjacentHTML('beforeend', entryHtml);
        });
    }
}

// Populate the auxiliary Enzyme forms, triggered by DOMContentLoaded event
function populateAuxEnzymeForm(container, data) {
    if (data.enzyme?.auxiliaryEnzymes) {
        data.enzyme.auxiliaryEnzymes.forEach(entry => {
            const index = container.children.length;
            const entryHtml = createHtmlAuxEnyzme(entry, index)
            container.insertAdjacentHTML('beforeend', entryHtml);
        });
    }
}

// Populate the Reaction forms, triggered by DOMContentLoaded event
function populateReactionForm(container, data, form_vals) {
    if (data.reactions) {
        data.reactions.forEach(entry => {
            const index = container.children.length;
            createHtmlReaction(entry, index, form_vals)
        });
    }
}

// Triggered by button to add an empty reference form
function insertEnzymeRefForm(container) {
    const entryHtml = createHtmlRef("", "enzyme_ref", "enzyme_ref[]")
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Triggered by button to add an empty AuxEnzyme
function insertAuxEnzymeForm(container) {
    const index = container.children.length;
    const entryHtml = createHtmlAuxEnyzme({}, index);
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Triggered by button to add an empty ReactionForm
function insertReactionForm(container, form_vals) {
    const index = container.children.length;
    createHtmlReaction({}, index, form_vals);
}

// Triggered by button to add an empty ReactionReferenceForm
function insertReactionRefForm(index) {
    const container = document.getElementById(`reaction[${index}]ref-field`);
    const entryHtml = createHtmlRef("", "reaction_ref", `reaction[${index}]ref[]`);
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Triggered by button to add an empty KnownReactionForm
function insertKnownReactionForm(index) {
    const knownReactionElements = document.querySelectorAll('.knownreaction');
    const containerReaction = document.getElementById(`reaction[${index}]knownreaction-field`);
    const index_inner = containerReaction.children.length;
    addHtmlKnownReaction({}, index, index_inner);
}

// Creates an reference field, used by enzyme information and Reaction fields
function createHtmlRef(data, className, jsonID) {
    return `
        <div class="${className}">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="${jsonID}" pattern="^doi:10.+" id="${jsonID}" class="form-control form-control-sm" placeholder="n/a" value="${data}" required title="A Digital Object Identifier describing the reference publication. Must begin with 'doi:10.'">
                        <label for="${jsonID}" class="form-label">Reference DOI<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.${className}')">Remove</button>
                </div>
            </div>
        </div>
    `;
}

// Creates the auxilary enzyme field; index is used to keep track of the instance that is opened
function createHtmlAuxEnyzme(data = {}, index) {
    let displayIndex = index + 1;
    return `
        <div class="aux_enzyme">
            <div class="card card-body mb-2">
                <div class="row d-flex align-items-center m-2 g-2">
                    <div class="col">
                        <h6 class="mb-0">Auxiliary Enzyme #${displayIndex}</h6>
                    </div>
                    <div class="col-auto mx-auto">
                        <button type="button" class="btn btn-danger" onclick="removeField(this, '.aux_enzyme')">Remove</button>
                    </div>
                </div>
                <div class="row m-2 g-2">
                    <div class="col">
                         <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]name" id="auxenzyme[${index}]name" class="form-control form-control-sm" placeholder="n/a" value='${data?.name ?? ""}' required title="The commonly used name of the auxiliary enzyme (e.g. 'McjB')">
                            <label for="auxenzyme[${index}]name" class="form-label">Name<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]description" id="auxenzyme[${index}]description" class="form-control form-control-sm" placeholder="n/a" maxlength="50" value='${data?.description ?? ""}' required title="Auxiliary enzyme type description (e.g. 'Peptidase')">
                            <label for="auxenzyme[${index}]description" class="form-label">Type<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                        </div>
                    </div>
                </div>
                <div class="row m-2 g-2">
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]uniprot" id="auxenzyme[${index}]uniprot" class="form-control form-control-sm" placeholder="n/a" value='${data?.databaseIds?.uniprot ?? ""}' title="UniProtKB or UniParc ID of the enzyme">
                            <label for="auxenzyme[${index}]uniprot" class="form-label">UniProt ID<span class="mandatory">**</span> <i class="bi bi-question-circle"></i></label>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]genpept" id="auxenzyme[${index}]genpept" class="form-control form-control-sm" placeholder="n/a" value='${data?.databaseIds?.genpept ?? ""}' title="NCBI GenPept ID of the enzyme">
                            <label for="auxenzyme[${index}]genpept" class="form-label">GenPept ID<span class="mandatory">**</span> <i class="bi bi-question-circle"></i></label>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]wikidata" id="auxenzyme[${index}]wikidata" class="form-control form-control-sm" placeholder="n/a" value='${data?.databaseIds?.wikidata ?? ""}' title="Wikidata QID of the enzyme">
                            <label for="auxenzyme[${index}]wikidata" class="form-label">Wikidata ID <i class="bi bi-question-circle"></i></label>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
}

// Creates the reaction field; index is used to keep track of the instance that is opened
function createHtmlReaction(data = {}, index, form_vals) {
    const reactionEntry = document.createElement('div')
    reactionEntry.classList.add('reaction')
    let displayIndex = index + 1;
    reactionEntry.innerHTML = `
        <div class="card card-body">
            <div class="row mb-2 align-items-center">
                <div class="col">
                    <h5 class="fw-semibold lh-2 m-0">Reaction #${displayIndex}</h5>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.reaction')">Remove</button>
                </div>
            </div>
            <div class="row">
                <div class="col">
                    <div class="card card-body">
                        <div class="row m-2 g-2">
                            <div class="col">
                                <h6 class="fw-semibold lh-2">General Information</h6>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]smarts" name="reaction[${index}]smarts" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.reactionSMARTS ?? ""}' required title="Reaction SMARTS describing substrate specificity and the introduced changes. Include any (covalent) co-factors the substrate(s)/product(s) may be attached to. Use e.g. Ketcher to draw the reaction SMARTS.">
                                    <label for="reaction[${index}]smart" class="form-label">Reaction SMARTS<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                                </div>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]description" name="reaction[${index}]description" class="form-control form-control-sm"  placeholder="n/a" type="text" value='${data?.description ?? ""}' required title="Brief description of the reaction.">
                                    <label for="reaction[${index}]description" class="form-label">Description<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                                </div>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]rhea" name="reaction[${index}]rhea" class="form-control form-control-sm" type="text" placeholder="n/a" value='${data?.databaseIds?.rhea ?? ""}' title="RHEA Knowledgebase ID of reaction (e.g. '32647')">
                                    <label for="reaction[${index}]rhea" class="form-label">RHEA ID <i class="bi bi-question-circle"></i></label>
                                </div>
                            </div>
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]ec" name="reaction[${index}]ec" pattern="^[0-9]+\..+" class="form-control form-control-sm" placeholder="n/a" type="text" value='${data?.databaseIds?.ec ?? ""}' title="The EC (Enzyme Commission) Number of the reaction (e.g. '1.14.19.59')">
                                    <label for="reaction[${index}]ec" class="form-label">EC Number <i class="bi bi-question-circle"></i></label>
                                </div>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <div class="card card-body">
                                    <div class="row">
                                        <div class="col">
                                            <div class="form-text mb-2">Tailoring reaction terms<span class="mandatory">*</span> <i class="bi bi-question-circle" data-bs-toggle="tooltip" title="Controlled vocabulary of tailoring reactions. Select the most applicable term(s)."></i></div>
                                        </div>
                                    </div>
                                    <div class="row" id="reaction[${index}]tailoring-field"></div>
                                </div>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <div class="card card-body">
                                    <div class="row">
                                        <div class="col">
                                            <div class="form-text mb-2">Experimental evidence qualifiers<span class="mandatory">*</span> <i class="bi bi-question-circle" data-bs-toggle="tooltip" title="Controlled vocabulary of qualifiers describing the experimental evidence for the reaction."></i></div>
                                        </div>
                                    </div>
                                    <div class="row" id="reaction[${index}]evidencecode-field"></div>
                                </div>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <div id="reaction[${index}]ref-field"></div>
                            </div>
                        </div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <button type="button" class="btn btn-secondary" onclick="insertReactionRefForm('${index}')">Add Reaction Reference</button>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="col">
                    <div class="card card-body">
                        <div class="row m-2 g-2">
                            <div class="col">
                                <h6 class="fw-semibold lh-2">Known Reactions</h6>
                            </div>
                        </div>
                        <div id="reaction[${index}]knownreaction-field"></div>
                        <div class="row m-2 g-2">
                            <div class="col">
                                <button type="button" class="btn btn-secondary my-2" onclick="insertKnownReactionForm(${index})">Add Known Reaction</button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
    const container = document.getElementById('reactions');
    container.appendChild(reactionEntry);

    const containerTailoring = document.getElementById(`reaction[${index}]tailoring-field`);
    for (let tailoring of form_vals.tailoring) {
        if (data.tailoring && data.tailoring.includes(tailoring) ) {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col-3');
            entryHtml.classList.add('form-text');
            entryHtml.classList.add('text-break');
            entryHtml.innerHTML = addTickedCheckbox(`reaction[${index}]tailoring[]`, tailoring);
            containerTailoring.appendChild(entryHtml)
        } else {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col-3');
            entryHtml.classList.add('form-text');
            entryHtml.classList.add('text-break');
            entryHtml.innerHTML =  addUntickedCheckbox(`reaction[${index}]tailoring[]`, tailoring);
            containerTailoring.appendChild(entryHtml)
        }
    }

    const containerEvidenceCode = document.getElementById(`reaction[${index}]evidencecode-field`);
    for (let evidence of form_vals.evidence) {
        if (data.evidence?.evidenceCode && data.evidence.evidenceCode.includes(evidence) ) {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col-3');
            entryHtml.classList.add('form-text');
            entryHtml.classList.add('text-break');
            entryHtml.innerHTML = addTickedCheckbox(`reaction[${index}]evidencecode[]`, evidence);
            containerEvidenceCode.appendChild(entryHtml);
        } else {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col-3');
            entryHtml.classList.add('form-text');
            entryHtml.classList.add('text-break');
            entryHtml.innerHTML =  addUntickedCheckbox(`reaction[${index}]evidencecode[]`, evidence);
            containerEvidenceCode.appendChild(entryHtml);
        }
    };

    const containerEvidenceRef = document.getElementById(`reaction[${index}]ref-field`);
    if (data.evidence?.references) {
        data.evidence.references.forEach(ref => {
            const entryHtml = document.createElement('div');
            entryHtml.innerHTML =  createHtmlRef(ref, "reaction_ref", `reaction[${index}]ref[]`);
            containerEvidenceRef.appendChild(entryHtml);
        });
    };

    if (data.reactions) {
        data.reactions.forEach(data_reaction => {
            const containerKR = document.getElementById(`reaction[${index}]knownreaction-field`);
            const indexKnownReaction = containerKR.children.length;
            addHtmlKnownReaction(data_reaction, index, indexKnownReaction);
        });
    };
}

// Adds a ticked checkbox with the appropriate value
function addTickedCheckbox(container_id, qualifier) {
    return `
        <div class="form-check">
            <input class="form-check-input" type="checkbox" id="${qualifier}" name='${container_id}' value="${qualifier}" checked>
            <label class="form-check-label" for="${qualifier}">${qualifier}</label>
        </div>
    `;
}

// Adds an unticked checkbox with the appropriate value
function addUntickedCheckbox(container_id, qualifier) {
    return `
        <div class="form-check">
            <input class="form-check-input" type="checkbox" id="${qualifier}" name='${container_id}' value="${qualifier}">
            <label class="form-check-label" for="${qualifier}">${qualifier}</label>
        </div>
    `;
}

// Inserts a KnownReaction Form into an existing Reaction form; outer index tracks the instance of the Reaction, the inner index the instance of the Known reaction inside the Reaction
function addHtmlKnownReaction(data, index_outer, index_inner) {
    let displayInnerIndex = index_inner + 1;
    const knownReactionEntry = document.createElement('div');
    knownReactionEntry.classList.add('knownreaction');
    knownReactionEntry.innerHTML = `
        <div class="card card-body">
            <div class="row align-items-center m-2 g-2">
                <div class="col">
                    <h6 class="lh-2 mb-0">Known Reaction #${displayInnerIndex}</h6>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.knownreaction')">Remove</button>
                </div>
            </div>
            <div class="row m-2 g-2">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="reaction[${index_outer}]knownreaction[${index_inner}]substrate" id="reaction[${index_outer}]knownreaction[${index_inner}]substrate" class="form-control form-control-sm" placeholder="n/a" value='${data?.substrate ?? ""}' required title="The reaction substrate(s) SMILES string. Covalently attached co-factors must be specified. Multiple substrates/non-covalent co-factors can be specified using dot notation ('substrate1.substrate2'). Use e.g. Ketcher to draw the SMILES.">
                        <label for="reaction[${index_outer}]knownreaction[${index_inner}]substrate" class="form-label">Substrate SMILES<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                    </div>
                </div>
            </div>
            <div class="row m-2 g-2">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="reaction[${index_outer}]knownreaction[${index_inner}]description" id="reaction[${index_outer}]knownreaction[${index_inner}]description" class="form-control form-control-sm" placeholder="n/a" value='${data?.description ?? ""}' title="Description of the reaction example">
                        <label for="reaction[${index_outer}]knownreaction[${index_inner}]description" class="form-label">Description <i class="bi bi-question-circle"></i></label>
                    </div>
                </div>
            </div>
            <div class="row m-2 g-2">
                <div class="col-auto">
                    <div class="form-check-inline">
                        <input class="form-check-input" value="True" type="radio" name="reaction[${index_outer}]knownreaction[${index_inner}]intermediate" id="reaction[${index_outer}]knownreaction[${index_inner}]intermediate1">
                        <label class="form-check-label" for="reaction[${index_outer}]knownreaction[${index_inner}]intermediate1">True</label>
                    </div>
                    <div class="form-check-inline">
                        <input class="form-check-input" value="False" type="radio" name="reaction[${index_outer}]knownreaction[${index_inner}]intermediate" id="reaction[${index_outer}]knownreaction[${index_inner}]intermediate2">
                        <label class="form-check-label" for="reaction[${index_outer}]knownreaction[${index_inner}]intermediate2">False</label>
                    </div>
                </div>
                <div class="col">
                    <div class="form-text">Intermediate product?<span class="mandatory">*</span> <i class="bi bi-question-circle" data-bs-toggle="tooltip" title="Is/are the product/s intermediates or the final products of the biosynthetic pathway?"></i></div>
                </div>
            </div>
            <div class="row gx-2 mx-2 mb-0">
                <div class="col">
                    <div id="reaction[${index_outer}]knownreaction[${index_inner}]products[]"></div>
                </div>
            </div>
            <div class="row m-2 g-2">
                <div class="col">
                    <button type="button" class="btn btn-secondary my-2" onclick="insertProductForm('reaction[${index_outer}]knownreaction[${index_inner}]products[]')">Add Product</button>
                </div>
            </div>
            <div class="row gx-2 mx-2 m-0">
                <div class="col">
                    <div id="reaction[${index_outer}]knownreaction[${index_inner}]forbiddenproducts[]"></div>
                </div>
            </div>
            <div class="row m-2 g-2">
                <div class="col">
                    <button type="button" class="btn btn-secondary my-2" onclick="insertForbiddenProductForm('reaction[${index_outer}]knownreaction[${index_inner}]forbiddenproducts[]')">Add Forbidden Product</button>
                </div>
            </div>
        </div>
    `;
    const containerReaction = document.getElementById(`reaction[${index_outer}]knownreaction-field`);
    containerReaction.appendChild(knownReactionEntry);

    //    populates product fields based on input data
    if (data.products) {
        data.products.forEach(product => {
            insertProductForm(`reaction[${index_outer}]knownreaction[${index_inner}]products[]`, product);
        });
    };

    //    populates forbidden product fields based on input data
    if (data.forbidden_products) {
        data.forbidden_products.forEach(product => {
            insertForbiddenProductForm(`reaction[${index_outer}]knownreaction[${index_inner}]forbiddenproducts[]`, product);
        });
    };

    //    populates the "intermediate" checkbox based on input value
    if (data.isIntermediate) {
        let radioButton = document.getElementById(`reaction[${index_outer}]knownreaction[${index_inner}]intermediate1`);
        radioButton.checked = true;
    } else {
        let radioButton = document.getElementById(`reaction[${index_outer}]knownreaction[${index_inner}]intermediate2`);
        radioButton.checked = true;
    };

}

// Inserts a knownProduct form inside the Known Reactions field
function insertProductForm(container_id, data = "") {
    const reactionProductEntry = document.createElement('div');
    reactionProductEntry.classList.add('reactionproducts');
    reactionProductEntry.innerHTML = `
        <div class="row align-items-center mb-1">
            <div class="col">
                <div class="form-floating">
                    <input type="text" name="${container_id}" id="${container_id}" class="form-control form-control-sm" placeholder="n/a" value='${data}' required title="A single reaction product SMILES string. Covalently attached co-factors must be specified. Dot notation is not permitted and each product/non-covalent co-factor must be specified in a separate field. Use e.g. Ketcher to draw the SMILES.">
                    <label for="${container_id}" class="form-label">Product SMILES<span class="mandatory">*</span> <i class="bi bi-question-circle"></i></label>
                </div>
            </div>
            <div class="col-auto mx-auto">
                <button type="button" class="btn btn-danger" onclick="removeField(this, '.reactionproducts')">Remove</button>
            </div>
        </div>
    `;
    const container = document.getElementById(container_id);
    container.appendChild(reactionProductEntry);
}

// Inserts a ForbiddenProduct form inside the Known Reactions field
function insertForbiddenProductForm(container_id, data = "") {
    const forbiddenProductEntry = document.createElement('div');
    forbiddenProductEntry.classList.add('forbiddenproducts');
    forbiddenProductEntry.innerHTML = `
        <div class="row d-flex align-items-center mb-1">
            <div class="col">
                <div class="form-floating">
                    <input type="text" name="${container_id}" id="${container_id}" class="form-control form-control-sm" placeholder="n/a" value='${data}' required title="A single forbidden reaction product SMILES string ('must not' result from reaction). Covalently attached co-factors must be specified. Dot notation is not permitted and each product/non-covalent co-factor must be specified in a separate field. Use e.g. Ketcher to draw the SMILES.">
                    <label for="${container_id}" class="form-label">Forbidden Product SMILES <i class="bi bi-question-circle"></i></label>
                </div>
            </div>
            <div class="col-auto mx-auto">
                <button type="button" class="btn btn-danger" onclick="removeField(this, '.forbiddenproducts')">Remove</button>
            </div>
        </div>
    `;
    const container = document.getElementById(container_id);
    container.appendChild(forbiddenProductEntry);
}