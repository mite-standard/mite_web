document.addEventListener('keydown', function(event) {
    if (event.key === "ArrowLeft") {
        document.getElementById('return').click();
    }
    if (event.key === "ArrowRight") {
        document.getElementById('forward').click();
    }
});