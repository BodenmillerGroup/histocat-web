"""Add aquisition table

Revision ID: b182e0878258
Revises: 7bed832a38b0
Create Date: 2019-04-10 15:31:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB

# revision identifiers, used by Alembic.
revision = 'b182e0878258'
down_revision = '7bed832a38b0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table('acquisition',
                    sa.Column('id', sa.Integer(), primary_key=True, index=True),
                    sa.Column('experiment_id', sa.Integer(), index=True),
                    sa.Column('name', sa.String(), nullable=False, index=True),
                    sa.Column('description', sa.String()),
                    sa.Column('location', sa.String(4096)),
                    sa.Column('meta', JSONB()),
                    sa.Column('created_at', sa.DateTime(), nullable=False),
                    )
    op.create_foreign_key(
        'fk_acquisition_experiment',
        'acquisition', 'experiment',
        ['experiment_id'], ['id'],
    )


def downgrade():
    op.drop_constraint('fk_acquisition_experiment')
    op.drop_table('acquisition')
