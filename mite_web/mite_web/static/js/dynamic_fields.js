function setupDynamicFields(containerSelector, addButtonSelector, removeButtonSelector, groupSelector, fieldHtml) {
    document.querySelector(addButtonSelector).addEventListener('click', function() {
        var container = document.querySelector(containerSelector);
        container.insertAdjacentHTML('beforeend', fieldHtml);
    });

    document.querySelector(containerSelector).addEventListener('click', function(event) {
        if (event.target.matches(removeButtonSelector)) {
            event.target.closest(groupSelector).remove();
        }
    });
}