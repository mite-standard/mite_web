function removeEntry(button, class_ref) {
    const targetDiv = button.closest(class_ref);
    if (targetDiv) {
        targetDiv.remove();
    }
}

function addEnzymeRef() {
    const container = document.getElementById('enzyme-ref');
    const entryHtml = `
        <div class="enzyme_ref">
            <div class="row d-flex align-items-center">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="enzyme_ref[]" id="enzyme_ref[]" class="form-control" aria-describedby="Enzyme Reference DOI Help"  value="" required>
                        <label for="enzyme_ref[]" class="form-label">Enzyme Reference DOI</label>
                    </div>
                </div>
                <div class="col">
                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.enzyme_ref')">Remove</button>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function addAuxEnzyme() {
    const container = document.getElementById('aux-enzyme');
    const entryHtml = `
        <div class="aux_enzyme">
            <div class="row d-flex align-items-center">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="enzyme_ref[]" id="enzyme_ref[]" class="form-control" aria-describedby="Enzyme Reference DOI Help"  value="" required>
                        <label for="enzyme_ref[]" class="form-label">Enzyme Reference DOI</label>
                    </div>
                </div>
                <div class="col">
                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.aux_enzyme')">Remove</button>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}