const modal = document.querySelector("#figureModal");
const modalImage = document.querySelector("#modalImage");
const modalTitle = document.querySelector("#modalTitle");

document.querySelectorAll(".figure-button").forEach((button) => {
  button.addEventListener("click", () => {
    modalImage.src = button.dataset.figure;
    modalImage.alt = button.dataset.title;
    modalTitle.textContent = button.dataset.title;
    if (typeof modal.showModal === "function") {
      modal.showModal();
    }
  });
});

modal?.addEventListener("click", (event) => {
  if (event.target === modal) {
    modal.close();
  }
});
