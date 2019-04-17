"""Add slide table

Revision ID: a596fd31b7bd
Revises: 7bed832a38b0
Create Date: 2019-04-10 15:30:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = 'a596fd31b7bd'
down_revision = '7bed832a38b0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table('slide',
                    sa.Column('id', sa.Integer(), primary_key=True, index=True),
                    sa.Column('experiment_id', sa.Integer(), index=True),
                    sa.Column('name', sa.String(), nullable=False, index=True),
                    sa.Column('description', sa.String()),
                    sa.Column('location', sa.String(4096)),
                    sa.Column('meta', JSONB()),
                    sa.Column('created_at', sa.DateTime(), nullable=False),
                    )
    op.create_foreign_key(
        'fk_slide_experiment',
        'slide', 'experiment',
        ['experiment_id'], ['id'],
    )


def downgrade():
    op.drop_constraint('fk_slide_experiment')
    op.drop_table('slide')
