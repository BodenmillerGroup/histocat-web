"""Add submission table

Revision ID: 485aa67ca586
Revises: 96df6178591e
Create Date: 2019-06-18 12:03:38.580966

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = '485aa67ca586'
down_revision = '96df6178591e'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'submission',
        sa.Column('id', sa.BigInteger(), primary_key=True, index=True),
        sa.Column('program', sa.String(), nullable=False, index=True),
        sa.Column('experiment_id', sa.Integer(), index=True),
        sa.Column('user_id', sa.Integer(), index=True),
        sa.Column('top_task_id', sa.BigInteger(), index=True),
        sa.Column('created_at', sa.DateTime(), nullable=False),
    )
    op.create_foreign_key(
        'fk_submission_experiment',
        'submission', 'experiment',
        ['experiment_id'], ['id'],
    )
    op.create_foreign_key(
        'fk_submission_user',
        'submission', 'user',
        ['user_id'], ['id'],
    )


def downgrade():
    op.drop_constraint('fk_submission_experiment')
    op.drop_constraint('fk_submission_user')
    op.drop_table('submission')
