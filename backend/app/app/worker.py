from app.core.celery_app import celery_app


@celery_app.task(acks_late=True)
def test_celery(word: str):
    return f"test task return {word}"


@celery_app.task(acks_late=True)
def import_slide(experiment_id: int, uri: str):
    return f"import slide task return {uri}"
