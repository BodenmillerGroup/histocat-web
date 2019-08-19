<template>
  <v-dialog v-model="dialog" persistent max-width="600px">
    <template v-slot:activator="{ on }">
      <v-btn small v-on="on">
        Create Dataset
      </v-btn>
    </template>
    <v-card>
      <v-card-title>Create Dataset</v-card-title>
      <v-card-text>
        <v-form v-model="valid" ref="form" lazy-validation>
          <v-layout column>
            <v-flex>
              <v-text-field
                label="Name"
                v-model="name"
                v-validate="'required'"
                data-vv-name="name"
                :error-messages="errors.collect('name')"
                required
              ></v-text-field>
            </v-flex>
            <v-flex>
              <v-text-field
                label="Description"
                v-model="description"
              ></v-text-field>
            </v-flex>
          </v-layout>
        </v-form>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn @click="submit" :disabled="!valid">
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-dialog>
</template>

<script lang="ts">
  import { datasetModule } from '@/modules/datasets';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class CreateDatasetDialog extends Vue {
    readonly datasetContext = datasetModule.context(this.$store);

    dialog = false;

    valid = false;
    name: string = '';
    description: string = '';

    async mounted() {
      this.reset();
    }

    reset() {
      this.name = '';
      this.description = '';
      this.$validator.reset();
    }

    cancel() {
      this.dialog = false;
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        await this.datasetContext.actions.createDataset({ name: this.name, description: this.description });
        this.dialog = false;
      }
    }
  }
</script>
