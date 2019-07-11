"""Add slide table

Revision ID: a596fd31b7bd
Revises: 7bed832a38b0
Create Date: 2019-04-10 15:30:16.499555

"""
import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects.postgresql import JSONB, UUID

# revision identifiers, used by Alembic.

revision = 'a596fd31b7bd'
down_revision = '7bed832a38b0'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'slide',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('experiment_id', sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), index=True),
        sa.Column('metaname', sa.String(4096), index=True),
        sa.Column('original_id', sa.Integer(), index=True),
        sa.Column('original_metadata', sa.String()),
        sa.Column('meta', JSONB()),
        sa.Column('location', sa.String(4096)),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )
    op.create_unique_constraint(
        'uq_experiment_slide_metaname',
        'slide',
        ['experiment_id', 'metaname']
    )


def downgrade():
    op.drop_table('slide')
    op.drop_unique_constraint('uq_experiment_slide_metaname')

